//! BLASTX's tabular reports (outfmt 6 and 7, the default `std` fields with a tab as any field list
//! or delimiter, `report::tabular::TabularSpec::std`) and outfmt 10 and 20: the adapter of the shared
//! writer (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query of a query batch in order it
//! writes the query's warnings, outfmt 7's comment lines and one row per HSP of the query's final
//! hit list that has a Seq-align (`report::has_seqalign`), in the list's order; outfmt 20's
//! header line comes with the prolog (`write_prolog`, before the search) and outfmt 7's count of
//! queries with the epilog (`report::write_epilog`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:811-834
//! What `SetFields` reads of a BLASTX HSP is its Seq-align (a Std-seg of the query's nucleotides
//! and the subject's residues): the query's coordinates in the record (start above the end on a
//! minus frame) and its frame, the subject's coordinates, the scores, the alignment's lengths and
//! gaps, and the two rows (the query's codons translated, `report::alignment`), from which it
//! counts the identities and, with the formatter's matrix (`GetAsciiProteinMatrix` of the search's
//! matrix), the positives again (`report::tabular::TabularWriter::row`). The subject frame is 0
//! (`sframe`, and the second number of `frames`), `sstrand` is `N/A` (a protein subject) and
//! `qcovus` is `N/A` (NCBI computes it for blastn only).

use std::io::Write;

use anyhow::Result;

use super::args::ResolvedOptions;
use super::empty_rows::{
    empty_row_scores, reported_row_range, rows_empty, subject_mapped_rows, tabular_empty_row,
    AlignQueryIds, MappedRows, SubjectLoc, TabularEmptyRow,
};
use super::input::FastaRecord;
use super::report::{alignment, has_seqalign, query_endpoints, query_warning, ShownSubjects};
use super::runtime::BatchResults;
use crate::algorithm::common::tabular::FailedValues;
use crate::cli::NativeError;
use crate::report::pairwise::EmptyRowsRefusal;
use crate::report::serial::seqalign::{mapped_empty_query_nucleotide, ungapped_interval};
use crate::report::serial::seqalign::{
    row_range_exception, std_seg_lengths_not_accurate, AlignRow,
};
use crate::report::tabular::{
    coverage_percent, taxdb_warning_line, union_length, DeflineIds, QueryCoveragePerHsp,
    SubjectLookup, TabularField, TabularHsp, TabularKind, TabularQuery, TabularRows, TabularSpec,
    TabularSubject, TabularTaxonomy, TabularWriter,
};

/// outfmt 7's `# BLASTX 2.17.0+` line.
const PROGRAM_VERSION: &str = concat!("BLASTX ", crate::wording::report_version!());

/// outfmt 7's line of an RPSTBLASTN run (`ResolvedOptions::rpstblastn`).
const RPSTBLASTN_VERSION: &str = concat!("RPSTBLASTN ", crate::wording::report_version!());

/// The program line of the run's tabular reports (`args::ProgramName`).
fn program_version(program: super::args::ProgramName) -> &'static str {
    if program == super::args::ProgramName::RPSTBLASTN {
        RPSTBLASTN_VERSION
    } else {
        PROGRAM_VERSION
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-135
/// The warning that each formatter of `spec` writes when it is made (`x_CheckTaxDB`), after
/// the report's lines before it are flushed; nothing when `taxdb.bti` is found or no plural name
/// is among the fields.
fn write_taxdb_warning(
    program: &str,
    spec: &TabularSpec,
    taxdb_found: bool,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.checks_taxdb() && !taxdb_found {
        writer.flush()?;
        diagnostics.write_all(taxdb_warning_line(program).as_bytes())?;
        diagnostics.flush()?;
    }
    Ok(())
}

/// Writes NCBI's `PrintProlog` of a custom tabular report, before the search: for outfmt 20 the
/// warning of the formatter it makes and the line of the field names; nothing for the other
/// kinds.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
/// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254
/// BLASTX writes its prologs before the search, as NCBI (outfmt 0's too): a refusal of the first
/// query batch comes after them.
pub(super) fn write_prolog(
    program: super::args::ProgramName,
    spec: &TabularSpec,
    taxdb_found: bool,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.kind != TabularKind::CsvWithHeader {
        return Ok(());
    }
    write_taxdb_warning(program.lower, spec, taxdb_found, writer, diagnostics)?;
    TabularWriter::new(spec, program_version(program), None).prolog(writer)?;
    Ok(())
}

/// The warning of the formatter that outfmt 7's `PrintEpilog` makes for its last line
/// (`report::write_epilog` writes the line).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
pub(super) fn write_epilog_warning(
    program: super::args::ProgramName,
    spec: &TabularSpec,
    taxdb_found: bool,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.kind == TabularKind::Commented {
        write_taxdb_warning(program.lower, spec, taxdb_found, writer, diagnostics)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908
/// The fields that make `SetFields` read the alignment (`SetCounts`) whatever the engine's count
/// of identities; `nident` does only when that count is above 0.
const COUNTS_FIELDS: [TabularField; 16] = [
    TabularField::QueryStart,
    TabularField::QueryEnd,
    TabularField::SubjectStart,
    TabularField::SubjectEnd,
    TabularField::AlignmentLength,
    TabularField::Gaps,
    TabularField::GapOpenings,
    TabularField::QuerySeq,
    TabularField::SubjectSeq,
    TabularField::Positives,
    TabularField::Mismatches,
    TabularField::PercentPositives,
    TabularField::PercentIdentical,
    TabularField::QueryFrame,
    TabularField::Btop,
    TabularField::SubjectStrand,
];

/// The values of one subject that its rows show (`TabularSubject` without the coverages).
struct SubjectValues<'a> {
    id: std::borrow::Cow<'a, [u8]>,
    accession: &'a [u8],
    accession_version: &'a [u8],
    gi: i64,
    all_ids: std::borrow::Cow<'a, [DeflineIds]>,
    titles: Option<Vec<&'a [u8]>>,
    length: usize,
    taxonomy: &'a TabularTaxonomy,
    /// How the lookups of a subject fail (no id: `SetSubjectId`), `None` where they do not.
    failed: Option<FailedValues<'a>>,
}

impl SubjectValues<'_> {
    /// The subject of the shared writer, with the names of the taxonomy database (`None`: every
    /// lookup fails), `qcovs` of its HSPs (-1 for `N/A`) and its rank or index `oid`.
    fn subject<'s>(
        &'s self,
        names: Option<&'s crate::report::tabular::TabularNames>,
        qcovs: i32,
        oid: usize,
    ) -> TabularSubject<'s> {
        TabularSubject {
            id: &self.id,
            accession: if self.accession.is_empty() {
                &self.id
            } else {
                self.accession
            },
            accession_version: if self.accession_version.is_empty() {
                &self.id
            } else {
                self.accession_version
            },
            gi: self.gi,
            all_ids: &self.all_ids,
            length: self.length,
            titles: self.titles.as_deref(),
            taxonomy: self.taxonomy,
            names,
            nucleotide: false,
            qcovs,
            qcovus: -1,
            lookup: self.lookup(),
            index: oid,
        }
    }

    /// The lookup of the subject's ids (`report::tabular::SubjectLookup`).
    fn lookup(&self) -> SubjectLookup<'_> {
        match &self.failed {
            Some(failed) => failed.lookup(&self.all_ids),
            None => SubjectLookup::Found,
        }
    }
}

/// The values of the subject `subject` (rank or index `oid`) of a report of `shown`; `slen` is the
/// length of the record, also with `-subject_loc`. The `sall*` of a record are made only for a
/// report that shows them (`all_deflines`, `TabularNeeds::all_deflines`).
fn subject_values<'a>(
    shown: &'a ShownSubjects<'_>,
    subject: &'a FastaRecord,
    oid: usize,
    no_taxonomy: &'a TabularTaxonomy,
    all_deflines: bool,
) -> SubjectValues<'a> {
    match shown {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:560-565
        // A FASTA subject has no definition lines: its `sall*` are its own ids (the first
        // word of its title), its title and taxonomy are empty (`N/A`).
        ShownSubjects::Fasta { parsed, .. }
            if parsed.tabular.get(oid).is_some_and(Option::is_some) =>
        {
            let ids = parsed.tabular[oid].as_ref().expect("the guard");
            // No ids: `SetSubjectId`'s `GetSeqIdList` throws (`RecordTabularIds::ids`); a
            // record's definition line is not looked up.
            let failed = parsed.subjects[oid]
                .as_deref()
                .and_then(|record| FailedValues::of_parsed(record, ids.ids.is_none(), false));
            let [id, accession, accession_version] = match &ids.ids {
                Some([id, accession, version]) => [id.as_str(), accession, version],
                None => [""; 3],
            };
            SubjectValues {
                id: std::borrow::Cow::Borrowed(id.as_bytes()),
                accession: accession.as_bytes(),
                accession_version: accession_version.as_bytes(),
                gi: parsed.subjects[oid]
                    .as_deref()
                    .map_or(0, |subject| subject.tabular_gi()),
                all_ids: std::borrow::Cow::Borrowed(&ids.all_ids),
                titles: None,
                length: subject.placement().length,
                taxonomy: no_taxonomy,
                failed,
            }
        }
        // A record without parsed ids: its `sseqid` is the first word of its title, or of the
        // title that `GenerateDefline` makes of it for a `Subject_` word
        // (`FastaRecord::subject_tabular_id`), none where that throws; its `sall*` the first word of
        // its title (`s_ReplaceLocalId`, `DeflineIds::of_record`).
        ShownSubjects::Fasta { .. } => {
            let tabular = subject.subject_tabular_id();
            let failed_id = tabular.is_none();
            let id = std::borrow::Cow::Owned(tabular.unwrap_or_default());
            SubjectValues {
                accession: &[],
                accession_version: &[],
                gi: 0,
                // Only where the word is not the id, which the writer prints for none.
                all_ids: if all_deflines && subject.shown_id() != &*id {
                    std::borrow::Cow::Owned(DeflineIds::of_record(subject.shown_id()))
                } else {
                    std::borrow::Cow::Borrowed(&[])
                },
                titles: None,
                length: subject.placement().length,
                taxonomy: no_taxonomy,
                failed: failed_id.then(|| {
                    let local = std::borrow::Cow::Borrowed(subject.internal_id.as_bytes());
                    FailedValues {
                        align_id: [local.clone(), local],
                        subject_id: true,
                        all_ids: false,
                        thrown: None,
                    }
                }),
                id,
            }
        }
        // A database sequence: the default report's ids (`ReportedSubjects::tabular_ids`), with
        // parsed Seq-ids the ids and titles of each definition line, else its one stored title.
        ShownSubjects::Database { reported, .. } => {
            // A subject whose load throws where the query's batch has its prefetch off
            // (`report::db_subject::HitLoads::thrown`): NCBI reads none of its own values.
            if let (Some(exception), Some(subject)) =
                (reported.loads().thrown(oid), reported.parsed_subject(oid))
            {
                return SubjectValues {
                    id: std::borrow::Cow::Borrowed(&[]),
                    accession: &[],
                    accession_version: &[],
                    gi: 0,
                    all_ids: std::borrow::Cow::Borrowed(&[]),
                    titles: None,
                    length: 0,
                    taxonomy: no_taxonomy,
                    failed: Some(FailedValues::thrown(subject, exception)),
                };
            }
            // No ids where the lookup fails (`ReportedSubjects::tabular_ids`).
            let found = reported.tabular_ids(oid);
            let [id, accession, accession_version] = found.unwrap_or([""; 3]);
            let parsed = reported.parsed_subject(oid);
            let (all_ids, all_failed) = reported.tabular_all_ids(oid);
            SubjectValues {
                id: std::borrow::Cow::Borrowed(id.as_bytes()),
                accession: accession.as_bytes(),
                accession_version: accession_version.as_bytes(),
                gi: parsed.map_or(0, |subject| subject.tabular_gi()),
                all_ids: std::borrow::Cow::Borrowed(all_ids),
                titles: Some(match parsed {
                    Some(subject) => subject
                        .deflines()
                        .iter()
                        .map(|defline| defline.title.as_deref().unwrap_or_default())
                        .collect(),
                    None => vec![reported
                        .stored_titles()
                        .get(oid)
                        .map_or(&[][..], |title| title)],
                }),
                length: subject.placement().length,
                taxonomy: reported.tabular_taxonomy(oid).unwrap_or(no_taxonomy),
                // No ids: `GetSeqIdList` throws, for the general id of a sequence without
                // parsed Seq-ids (`TabularSubjectId::Unknown`) in `SetSubjectId` and
                // `x_SetSubjectIds`; for one with parsed Seq-ids where `ParsedSubject` says.
                failed: match parsed {
                    Some(parsed) => FailedValues::of_parsed(parsed, found.is_none(), all_failed),
                    None => found.is_none().then(|| {
                        let [full, accession] =
                            crate::algorithm::common::database::failed_general_lookup(
                                reported.oid(oid),
                            );
                        FailedValues {
                            align_id: [
                                std::borrow::Cow::Owned(full),
                                std::borrow::Cow::Owned(accession),
                            ],
                            subject_id: true,
                            all_ids: true,
                            thrown: None,
                        }
                    }),
                },
            }
        }
    }
}

/// Writes the custom tabular report `spec` of the query batches `batches` (their final hit
/// lists sorted by e-value, `report::render_queries`) to `writer`: for each query its warning
/// (`report::query_warning`) and the warning of its formatter on `stderr`, outfmt 7's lines and
/// the rows. A database's refusals of the fields come before (`ReportedSubjects::read_tabular`,
/// `native::run`); each row's texts of the taxonomy database's lookups follow the rows before
/// it (`report::tabular::TabularNames::row_messages`). An error writing a query's report is the
/// stream failure of NCBI's formatter (`ostream::ReportWriteFailure`), at the flush of each
/// query's `CBlastTabularInfo` where its writes were held.
#[allow(clippy::too_many_arguments)]
pub(super) fn write_blastx_tabular(
    writer: &mut crate::cli::NcbiOstream<'_>,
    stderr: &mut dyn Write,
    options: &ResolvedOptions,
    spec: &TabularSpec,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
    mut probe: Option<&mut crate::api::local_blast::FormatProbe<'_>>,
    first_hsp: crate::api::local_blast::HspIndex,
) -> Result<()> {
    if spec.kind == crate::report::tabular::TabularKind::TaxReport {
        return write_blastx_tax_report(writer, stderr, options, queries, subjects, shown, batches);
    }
    let needs = spec.needs();
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:183
    // The search's matrix (`-matrix`, item SC6 of line NP). RPSTBLASTN's options carry no matrix
    // name (the `.aux` matrix goes to the search's options only), so NCBI's formatter has no
    // matrix and counts the identities as the positives, as rpsblast's
    // (`blastp::tabular::write_rpsblast_tabular`).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:238
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3552-3556
    let table = TabularWriter::new(
        spec,
        program_version(options.program()),
        options.rpstblastn.is_none().then_some(options.scoring),
    );
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1311
    let database: Vec<u8> = match shown {
        ShownSubjects::Fasta { path, .. } => crate::report::pairwise::subject_set_name(path),
        ShownSubjects::Database { name, .. } => name.as_bytes().to_vec(),
    };
    let names = match shown {
        ShownSubjects::Database { reported, .. } => reported.tabular_names(),
        ShownSubjects::Fasta { .. } => None,
    };
    let counts_by_fields = COUNTS_FIELDS.iter().any(|field| spec.has(*field));
    // Which rows read their subject (`TabularSpec::subject_fetch`).
    let fetch = spec.subject_fetch(crate::blastinput::app::Program::Blastx);
    let no_taxonomy = TabularTaxonomy::default();
    let loc = SubjectLoc::of(shown);
    let mut query_ranges: Vec<(u64, u64)> = Vec::new();
    let mut next_hsp = first_hsp;
    // The bytes of one comment block or one row, written to the report stream as soon as they
    // are made (as NCBI's insertions reach its stream): nothing of a query's report is held.
    let mut buf: Vec<u8> = Vec::new();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
            stderr.write_all(&query_warning(options, query, batch.search_skipped))?;
            stderr.write_all(shown.query_fetches(batch.query_ordinal + qidx))?;
            // The texts of the subjects that the query's rows read one by one, written where they
            // read each (`ReadCursor`, the query's only).
            let mut reads = crate::report::query_warnings::ReadCursor::new(std::slice::from_ref(
                shown.query_reads(batch.query_ordinal + qidx),
            ));
            // NCBI makes a formatter for each query's report, and it warns when the field list
            // asks for a plural name of the taxonomy database that it cannot find
            // (tabular.cpp:139-158, `x_CheckTaxDB`).
            write_taxdb_warning(
                options.program().lower,
                spec,
                options.taxdb_found,
                &mut *writer,
                &mut *stderr,
            )?;
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
            // A stream failure unwinds through the query's formatter.
            (|| -> Result<()> {
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1277-1278
                // A query of a batch that NCBI did not search has no Seq-align set: no count.
                // (RPSTBLASTN: also the HSPs of the profiles that left the merged lists past the
                // formatter's cut, `BatchResults::unlisted_hsps`.)
                let count: usize = lists
                    .iter()
                    .map(|list| {
                        list.hsps
                            .iter()
                            .filter(|hsp| has_seqalign(hsp, options.gapped))
                            .count()
                    })
                    .sum::<usize>()
                    + batch.unlisted_hsps(qidx);
                table.query_header(
                    &mut buf,
                    &query.query_line(),
                    &database,
                    // A query with an alignment set (`BatchResults::has_alignment_set`).
                    batch.has_alignment_set(options).then_some(count),
                )?;
                writer.write_all(&buf)?;
                buf.clear();
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:787-790
                // The searched letters: the record, or its `-query_loc` interval (`qcovhsp` below);
                // `m_QueryRange` is the range as given (`query_length`).
                let searched_length = query.sequence.len();
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:880-887
                // `qlen` and the `qcovs` of the subjects: the length of the `-query_loc` range as
                // given (even where it ends past the record), else the query's.
                let query_length = options
                    .query_range
                    .map_or(searched_length, |range| range.length());
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1497-1499
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1034-1043
                // The length against which `hsp_percent_coverage` counts an HSP's residues in its
                // frame (`qcovhsp`): a third of the searched letters in a gapped search, the
                // searched letters in an ungapped one (blast_seqalign.cpp:1382-1383).
                let coverage_length = if options.gapped {
                    searched_length as i32 / 3
                } else {
                    searched_length as i32
                };
                let tabular_query =
                    TabularQuery::of_ids(query.believed(), query.shown_id(), query_length);
                // Several volumes: the merged lists cut once per query (`BatchResults::shown_lists`).
                let listed = batch.shown_lists(qidx, options);
                // The query's formatter (`TabularRows`).
                let mut table_rows = TabularRows::default();
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:816-824
                // `qcovs` reads the query range of every HSP of the query before its first row
                // (`GetMasterCoverage`, cited below): a query row whose range it cannot take
                // throws there, after outfmt 7's lines (`empty_rows::reported_row_range`: the
                // empty row of an HSP of length 0 with `-db`, with `-subject` where the mapper
                // leaves it or its mapped query row wraps (`MappedRows::QueryRowWraps`), and a
                // leading gap left in a reported gapped edit script for `-db` and `-subject`
                // alike).
                let unranged = needs
                    .qcovs
                    .then(|| {
                        listed.iter().find_map(|list| {
                            let list_loc = loc.of_list(list, options.gapped);
                            let empty_query_rewritten = list_loc == SubjectLoc::Interval
                                && !matches!(
                                    subject_mapped_rows(
                                        list_loc,
                                        list,
                                        options.gapped,
                                        query,
                                        &subjects[list.oid]
                                    ),
                                    Some(MappedRows::QueryRowWraps)
                                );
                            list.hsps.iter().find_map(|hsp| {
                                has_seqalign(hsp, options.gapped)
                                    .then(|| {
                                        reported_row_range(
                                            hsp,
                                            AlignRow::Query,
                                            empty_query_rewritten,
                                        )
                                    })
                                    .flatten()
                            })
                        })
                    })
                    .flatten();
                if let Some(error) = unranged {
                    writer.flush()?;
                    return Err(NativeError::exception(&row_range_exception(error)).into());
                }
                for list in listed {
                    let shown_hsps = || {
                        list.hsps
                            .iter()
                            .filter(|hsp| has_seqalign(hsp, options.gapped))
                    };
                    let shown_count = shown_hsps().count();
                    if shown_count == 0 {
                        continue;
                    }
                    let subject = &subjects[list.oid];
                    let values =
                        subject_values(shown, subject, list.oid, &no_taxonomy, needs.all_deflines);
                    // The subject's Seq-loc as its Seq-align shows it (`SubjectLoc::of_list`).
                    let loc = loc.of_list(list, options.gapped);
                    // The subject's `qcovhsp` rule: its rows only divide.
                    let qcovhsp = QueryCoveragePerHsp::of_subject(
                        !options.gapped,
                        shown_count,
                        coverage_length,
                        query_length as i32,
                    );
                    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:629-640
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2969-2981
                    // The HSPs of one subject are one group of `InsertSubjectScores` (a hit list
                    // holds a subject once): the union of their query intervals in the record's
                    // nucleotides, in percent of the query's length.
                    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1541-1547
                    // A `-subject` subject whose only shown HSP has length 0: the mapper's rows
                    // without scores (`empty_rows::MappedRows`, which `render_queries_with` let
                    // through).
                    let mapped = (shown_count == 1)
                        .then(|| subject_mapped_rows(loc, list, options.gapped, query, subject))
                        .flatten();
                    let qcovs = if needs.qcovs {
                        query_ranges.clear();
                        query_ranges.extend(shown_hsps().map(|hsp| {
                            let raw = &hsp.hsp;
                            if let Some(MappedRows::Letters(mapped)) = &mapped {
                                let (start, end) = (mapped.q_start as u64, mapped.q_end as u64);
                                return (start.min(end) - 1, start.max(end) - 1);
                            }
                            // The one nucleotide of `-subject`'s mapped empty row
                            // (`seqalign::mapped_empty_query_nucleotide`); `-db` stopped above.
                            if rows_empty(hsp) {
                                let at = mapped_empty_query_nucleotide(
                                    ungapped_interval(
                                        raw.q_start,
                                        raw.q_end,
                                        raw.frame,
                                        searched_length as i32,
                                    )
                                    .0,
                                    raw.frame,
                                ) as u64;
                                return (at, at);
                            }
                            let (start, end) = query_endpoints(
                                searched_length,
                                raw.frame,
                                raw.q_start,
                                hsp.seqalign_ends().0,
                            );
                            (start.min(end) as u64 - 1, start.max(end) as u64 - 1)
                        }));
                        coverage_percent(union_length(&mut query_ranges), query_length)
                    } else {
                        -1
                    };
                    let tabular_subject = values.subject(names, qcovs, list.oid);
                    // A subject whose lookups fail shows what the last found subject set
                    // (`TabularRows::last_found`), made again from its index.
                    let kept_values = match table_rows.last_found() {
                        Some(kept) if table.fails(&tabular_subject) => Some(subject_values(
                            shown,
                            &subjects[kept],
                            kept,
                            &no_taxonomy,
                            needs.all_deflines,
                        )),
                        _ => None,
                    };
                    let kept = kept_values
                        .as_ref()
                        .zip(table_rows.last_found())
                        .map(|(values, kept)| values.subject(names, -1, kept));
                    for hsp in shown_hsps() {
                        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
                        // A row that reads the subject (`TabularSpec::subject_fetch`, by the
                        // engine's identities for `nident` alone) loads it first where the
                        // query's batch did not prefetch it: the texts of its taxonomy lookups
                        // follow the rows before it, before those of the row's own lookups.
                        // A row that reads a subject whose load throws outside the subject block
                        // ends the run (`TabularWriter::ends_run`).
                        if let Some(reads) = reads.as_mut() {
                            if fetch.reads(hsp.identity) {
                                reads.subject_read(
                                    0,
                                    list.oid as u32,
                                    &mut *writer,
                                    &mut *stderr,
                                )?;
                            }
                        }
                        if let Some(end) = table.ends_run(
                            &tabular_subject,
                            crate::blastinput::app::Program::Blastx,
                            hsp.identity,
                        ) {
                            return Err(crate::report::tabular::run_end(writer, &buf, end, || {
                                if let Some(reads) = reads.as_mut() {
                                    reads.end_run();
                                }
                            }));
                        }
                        // The texts of the row's lookups of the taxonomy database
                        // (`x_SetTaxInfoAll`, `x_SetTaxInfo`) follow the rows before it.
                        if let (Some(names), true) = (names, table.sets_subject(&tabular_subject)) {
                            let text = names.row_messages(spec, values.taxonomy);
                            if !text.is_empty() {
                                writer.flush()?;
                                stderr.write_all(&text)?;
                                stderr.flush()?;
                            }
                        }
                        let raw = &hsp.hsp;
                        // An HSP of length 0 (`empty_rows::tabular_empty_row`).
                        let empty = rows_empty(hsp);
                        let form = if empty {
                            tabular_empty_row(
                                loc,
                                spec.opens_counts(&hsp.identity),
                                shown_count,
                                spec.reads_query(),
                                matches!(mapped, Some(MappedRows::QueryRowWraps)),
                            )
                        } else {
                            TabularEmptyRow::Row
                        };
                        let (evalue, bit_score, score) = if empty {
                            empty_row_scores(loc, shown_count, hsp)
                        } else {
                            (hsp.evalue, hsp.bit_score, raw.score)
                        };
                        if form == TabularEmptyRow::Throws {
                            writer.flush()?;
                            return Err(
                                NativeError::exception(&std_seg_lengths_not_accurate()).into()
                            );
                        }
                        let a;
                        let catch_ids;
                        let (row, row_query, row_subject) = if form == TabularEmptyRow::Mapped {
                            // `render_queries_with` refused the subjects whose alignment LOSAT
                            // does not map (`empty_rows::refuse_mapped_without_scores`).
                            let Some(MappedRows::Letters(mapped)) = &mapped else {
                                return Err(EmptyRowsRefusal(
                                    super::empty_rows::UNPORTED_SUBJECT_ALL_EMPTY.to_string(),
                                )
                                .into());
                            };
                            // No scores (`empty_row_scores`), and the identities of `nident`
                            // outside the alignment block `x_ResetFields`' 0 (`GetAlnScores`
                            // finds no `num_ident`); inside it every count is of the rows' letters.
                            let row = TabularHsp {
                                q_start: mapped.q_start,
                                q_end: mapped.q_end,
                                s_start: mapped.s_start,
                                s_end: mapped.s_end,
                                evalue,
                                bit_score,
                                score,
                                length: mapped.length(),
                                gap_opens: 0,
                                gaps: 0,
                                num_ident: 0,
                                q_frame: mapped.q_frame,
                                s_frame: 0,
                                rows: Some(mapped.rows()).filter(|_| counts_by_fields),
                                // The mapped row's `-` past the query's end is the gap character
                                // inside its segment, not an end gap (alnvec.cpp:279-297,
                                // `TabularHsp::ops`): no script.
                                ops: &[],
                                qcovhsp: qcovhsp.of_hsp(0, mapped.length() as i32),
                            };
                            (row, tabular_query, tabular_subject)
                        } else if form == TabularEmptyRow::Leaves {
                            // `slen` is 0 without a subject field that reads its Bioseq
                            // (`TabularSpec::reads_subject_bioseq`).
                            let row = TabularHsp {
                                q_start: 0,
                                q_end: 0,
                                s_start: 0,
                                s_end: 0,
                                evalue,
                                bit_score,
                                score,
                                length: 0,
                                gap_opens: 0,
                                gaps: 0,
                                num_ident: 0,
                                q_frame: 0,
                                s_frame: 0,
                                rows: None,
                                ops: &[],
                                qcovhsp: -1,
                            };
                            catch_ids = AlignQueryIds::of(query);
                            let catch_query = TabularQuery {
                                id: &catch_ids.ids[0],
                                accession: &catch_ids.ids[1],
                                accession_version: &catch_ids.ids[2],
                                gi: catch_ids.gi,
                                length: query_length,
                            };
                            let catch_subject = TabularSubject {
                                length: if spec.reads_subject_bioseq() {
                                    values.length
                                } else {
                                    0
                                },
                                ..tabular_subject
                            };
                            (row, catch_query, catch_subject)
                        } else {
                            // The writer counts the rows' identities and positives
                            // (`TabularWriter::row`); `alignment` makes no counts.
                            a = alignment(hsp, subject, query, &code)?;
                            // NCBI reads the alignment's rows when a field reads its counts
                            // (`COUNTS_FIELDS`, or `nident` with the engine's identities above 0;
                            // BLASTX fetches the sequences, so `nident`, `mismatch` and `pident`
                            // count the rows' identities); otherwise `nident` is the engine's count.
                            let rows = (counts_by_fields || hsp.identity > 0)
                                .then_some((a.qseq.as_bytes(), a.sseq.as_bytes()));
                            let row = TabularHsp {
                                q_start: i64::from(a.qstart),
                                q_end: i64::from(a.qend),
                                s_start: i64::from(a.sstart),
                                s_end: i64::from(a.send),
                                evalue,
                                bit_score,
                                score,
                                length: a.qseq.len(),
                                gap_opens: a.gap_opens,
                                gaps: a.gaps,
                                num_ident: hsp.identity,
                                q_frame: i32::from(a.frame),
                                s_frame: 0,
                                rows,
                                ops: &hsp.edit_script,
                                qcovhsp: qcovhsp.of_hsp(
                                    // The HSP's own ends (`Blast_HSPGetQueryCoverage` of the
                                    // BlastHSP, blast_seqalign.cpp:1223-1224), not the Seq-align's,
                                    // with the letters that `-subject`'s mapper left out.
                                    raw.q_end - raw.q_start + hsp.query_left_out as i32,
                                    // The Seq-align's nucleotides; `-subject`'s mapped empty row
                                    // has one.
                                    if empty {
                                        1
                                    } else {
                                        (a.qstart - a.qend).abs() + 1
                                    },
                                ),
                            };
                            (row, tabular_query, tabular_subject)
                        };
                        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1102-1108
                        // One printed row is one HSP; the probe marks it without changing it.
                        if let Some(probe) = probe.as_deref_mut() {
                            writer.flush()?;
                            probe.begin(next_hsp);
                        }
                        table.row(
                            &mut buf,
                            &mut table_rows,
                            &row_query,
                            &row_subject,
                            kept.as_ref(),
                            &row,
                        );
                        writer.write_all(&buf)?;
                        buf.clear();
                        if let Some(probe) = probe.as_deref_mut() {
                            writer.flush()?;
                            probe.end(next_hsp);
                        }
                        next_hsp += 1;
                    }
                }
                // A text whose read the rows did not place comes after them.
                if let Some(reads) = reads.as_mut() {
                    reads.before_query(1, &mut *writer, &mut *stderr)?;
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
                crate::report::ostream::ReportWrite::ncbi_flush(writer)?;
                Ok(())
            })()?;
        }
    }
    Ok(())
}

/// Writes the organism report (outfmt 18, `report::tax_report`) of the query batches `batches`
/// to `writer`: for each query its warning (`report::query_warning`) on `stderr`, then the
/// caption, the `Query=` block (the query record's length, also with `-query_loc`, as outfmt 0's)
/// and the table of the HSPs that the tabular rows show (`write_blastx_tabular`), the texts of
/// the subjects that the report reads one by one where it reads each (`ShownSubjects::query_reads`,
/// `TaxReport::write_query`). The report is made before the batches' first byte
/// (`tax_report::report_batch`), so that LOSAT's refusals come before it. NCBI writes it through
/// `m_Outfile` (`CBlastFormat`), not a `CBlastTabularInfo`: a failing write is outfmt 0's
/// "BLAST failed to write output" (`native::report_error` of the `io::Error`), not the tabular
/// formats' abort (probe `T/t2/fix1` `f1.d5.blastx.*`).
#[allow(clippy::too_many_arguments)]
fn write_blastx_tax_report(
    writer: &mut crate::cli::NcbiOstream<'_>,
    stderr: &mut dyn Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
) -> Result<()> {
    use crate::report::tax_report::{report_batch, TaxHsp, TaxPlain, TaxQuery};
    let program = options.program().upper;
    let titles: Vec<Vec<u8>> = batches
        .iter()
        .flat_map(|batch| {
            (0..batch.queries.len()).map(|qidx| {
                queries[batch.query_ordinal + qidx]
                    .query_line()
                    .into_owned()
            })
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:893-899
    // The HSPs of every list of each query: NCBI does not cut the organism report's Seq-aligns
    // to the hit list size (no `PruneSeqalign` or `TrimSeqAlign`, as for the archive), so
    // RPSTBLASTN's merged lists of several volumes are read whole (`BatchResults::shown_lists`
    // is the cut of the other reports).
    let tax_queries: Vec<TaxQuery<'_>> = batches
        .iter()
        .flat_map(|batch| {
            batch.queries.iter().enumerate().map(move |(qidx, lists)| {
                (
                    batch.query_ordinal + qidx,
                    lists.as_slice(),
                    batch.has_alignment_set(options),
                )
            })
        })
        .zip(&titles)
        .map(|((ordinal, lists, alignment_set), title)| TaxQuery {
            title,
            length: queries[ordinal].placement().length,
            // A query without an alignment set (`BatchResults::has_alignment_set`: of a batch
            // that NCBI did not search, but for the merge of several volumes) has no
            // Seq-align-set.
            hsps: alignment_set.then(|| {
                lists
                    .iter()
                    .flat_map(|list| {
                        list.hsps
                            .iter()
                            .filter(|hsp| has_seqalign(hsp, options.gapped))
                            .map(move |hsp| TaxHsp {
                                subject: list.oid,
                                evalue: hsp.evalue,
                                bit_score: hsp.bit_score,
                            })
                    })
                    .collect()
            }),
        })
        .collect();
    let parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>] = match shown {
        ShownSubjects::Fasta { parsed, .. } => &parsed.subjects,
        ShownSubjects::Database { reported, .. } => reported.parsed(),
    };
    let plain = |s_idx: usize| match shown {
        ShownSubjects::Fasta { .. } => TaxPlain::Record {
            local_id: subjects[s_idx].internal_id.as_bytes(),
            title: &subjects[s_idx].title,
        },
        ShownSubjects::Database { reported, .. } => TaxPlain::Stored {
            oid: reported.oid(s_idx),
            title: reported
                .stored_titles()
                .get(s_idx)
                .map_or(&[][..], |title| title),
            source: crate::report::pairwise::stored_source(reported.sources(), s_idx as u32),
            taxid: reported
                .tabular_taxonomy(s_idx)
                .map_or_else(|| reported.taxid(s_idx), |taxonomy| taxonomy.taxid),
        },
    };
    // A `-subject` search's report looks names up in the run's database.
    let subject_names = options
        .subject_taxdb
        .as_ref()
        .map(|taxdb| taxdb.report_names(program));
    let names = match shown {
        ShownSubjects::Database { reported, .. } => reported.tabular_names(),
        ShownSubjects::Fasta { .. } => subject_names.as_ref(),
    };
    let loads = match shown {
        ShownSubjects::Database { reported, .. } => Some(reported.loads()),
        ShownSubjects::Fasta { .. } => None,
    };
    let report = report_batch(
        &tax_queries,
        parsed,
        plain,
        true,
        loads,
        names,
        program,
        options.is_html,
    )?;
    let mut q_idx = 0;
    for batch in batches {
        for qidx in 0..batch.queries.len() {
            let query = &queries[batch.query_ordinal + qidx];
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
            // `cerr` is tied to `cout`: the report before the texts reaches the output first.
            let warning = query_warning(options, query, batch.search_skipped);
            if !warning.is_empty() {
                writer.flush()?;
            }
            stderr.write_all(&warning)?;
            // The texts of the subjects that the query's report reads one by one, written where
            // it reads each (`ReadCursor`, the query's only).
            let mut reads = crate::report::query_warnings::ReadCursor::new(std::slice::from_ref(
                shown.query_reads(batch.query_ordinal + qidx),
            ));
            let tax_query = &tax_queries[q_idx];
            report.write_query(
                q_idx,
                &mut *writer,
                &mut *stderr,
                reads.as_mut().map(|cursor| (cursor, 0)),
                tax_query.title,
                tax_query.length,
            )?;
            // A text whose read the report did not place comes after it.
            if let Some(reads) = reads.as_mut() {
                reads.before_query(1, &mut *writer, &mut *stderr)?;
            }
            q_idx += 1;
        }
    }
    Ok(())
}
