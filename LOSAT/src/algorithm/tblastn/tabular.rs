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
//! What `SetFields` reads of an HSP is its Seq-align (`blast_xml` has the coordinates): the protein
//! query ascending, the subject in nucleotides (descending on a minus frame), the score, the
//! e-value (the sum e-value of a linked set) and the bit score, the alignment's lengths and gaps,
//! and the two rows of the alignment, which the TBLASTN formatter reads (no `SetNoFetch`) for
//! `qseq`, `sseq`, `btop`, `positive`, `ppos` and the identities: the protein query, unmasked, and
//! the subject translated with `-db_gencode` (the XML's rows without its masks,
//! `aligned_sequences_into`). Without those fields the identities are the default row's
//! (`report_num_ident`, counted from the same rows). `qcovus` is BLASTN's alone (`N/A`).

use std::io::Write;

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
use crate::cli::NativeError;
use crate::config::ScoringMatrix;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::query_warnings::QueryWarnings;
use crate::report::row_letters::QueryLetters;
use crate::report::serial::seqalign::{reported_std_seg_row_range, row_range_exception, AlignRow};
use crate::report::tabular::{
    coverage_percent, union_length, QueryCoveragePerHsp, TabularField, TabularHsp, TabularKind,
    TabularNames, TabularQuery, TabularRows, TabularSpec, TabularTaxonomy, TabularWriter,
};
use crate::utils::genetic_code::GeneticCode;

const PROGRAM: &str = "TBLASTN";
/// The `# TBLASTN 2.17.0+` line of outfmt 7 (`NStr::ToUpper` of `m_Program`, a space and `version`).
const PROGRAM_VERSION: &str = concat!("TBLASTN ", crate::wording::report_version!());
/// The `# PSITBLASTN 2.17.0+` line of a PSSM search (`-in_pssm`, `in_pssm::PROGRAM`).
const PSI_PROGRAM_VERSION: &str = concat!(
    super::in_pssm::psi_tblastn_name!(),
    " ",
    crate::wording::report_version!()
);

/// The program line of outfmt 7 (and 20's formatter): tblastn's, or PSI-TBLASTN's with
/// `-in_pssm` (`pssm`).
pub(super) fn program_version(pssm: bool) -> &'static str {
    if pssm {
        PSI_PROGRAM_VERSION
    } else {
        PROGRAM_VERSION
    }
}

/// Whether `nident` is the only field of the counts' block of `spec`, which reads the rows of some
/// HSPs (`write_tblastn_tabular`).
fn nident_alone(spec: &TabularSpec) -> bool {
    spec.has(TabularField::NumIdentical)
        && !spec.fields.iter().any(|field| COUNTS_BLOCK.contains(field))
}

/// Whether a report of `spec` can read the rows of an HSP (`write_tblastn_tabular`: the letters of
/// its subject): its fields read them, or `nident_alone`.
pub(super) fn reads_rows(spec: &TabularSpec) -> bool {
    spec.needs().rows || nident_alone(spec)
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
    /// An ungapped search (`-ungapped`), whose Seq-aligns are Std-segs (`qcovhsp`,
    /// `report::tabular::QueryCoveragePerHsp`).
    pub ungapped: bool,
    /// The `-html` links of the outfmt 0 to 4 reports (`blastinput::app::html_links`; their
    /// pairwise writer reads them here too, `stage_e_report::write_pairwise`); None without them.
    pub html: Option<std::sync::Arc<crate::report::html_ids::HtmlLinks>>,
    /// `-html` (`TblastnArgs::html`, NCBI's `m_IsHTML`): the HTML form of the outfmt 0 to 4 and
    /// 18 reports (the organism report takes no links).
    pub is_html: bool,
}

/// One custom tabular report of a run: `render`'s inputs that it reads.
pub(super) struct TblastnTabular<'a> {
    pub spec: &'a TabularSpec,
    pub inputs: &'a TabularInputs,
    /// The queries (their searched letters, which the rows read).
    pub queries: &'a [FastaRecord],
    /// The PSSM query's Bioseq letters where they are not its record's
    /// (`PsiSearch::bioseq_letters`), which the rows read (`row_letters::QueryLetters`); a run of
    /// such a query has no other.
    pub query_bioseq: Option<&'a [u8]>,
    /// The letters of the subjects by subject index (`-subject`) or rank (`-db`; their lengths
    /// alone where `reads_rows` is false).
    pub letters: SubjectLetters<'a>,
    /// The default report's tabular subject ids (`sseqid`, `sacc`, `saccver`), by subject index
    /// or rank (`stage_e_report::tabular_subject_ids`, `ReportedSubjects::tabular_subject_ids`).
    pub subject_ids: &'a [crate::algorithm::common::tabular::SubjectId],
    /// The records of a subject set (`-subject`; the organism report shows their local ids and
    /// titles).
    pub subject_records: &'a [FastaRecord],
    pub shown: &'a ShownSubjects<'a>,
    /// By query: whether NCBI did not search its batch (outfmt 7's lines without a count).
    pub query_batch_skipped: &'a [bool],
    pub ranges: &'a ReportRanges<'a>,
    /// `-db_gencode`: the code that translates the subject rows (`SetDbGeneticCode`).
    pub genetic_code: u8,
    /// The run's matrix (`m_ScoringMatrix`, the positives).
    pub matrix: ScoringMatrix,
    /// The program line of outfmt 7 (`program_version`).
    pub program_version: &'static str,
    /// A PSSM query of the molecule `other` (`PsiSearch::query_other_coding`): the first row
    /// whose fields read the query's letters throws (`report::tabular::query_sequence_end`).
    pub query_other_coding: Option<i32>,
}

/// The custom tabular report of a checked `-outfmt` value (`TabularSpec::of_app`, tblastn's
/// formatter, PSI-TBLASTN's too, gets the delimiter): a field list or a delimiter of outfmt 6 and
/// 7, and every outfmt 10 and 20; `None` for the default 6 and 7 (the `std` fields with a tab,
/// `stage_e_report::write_tabular`) and the other formats.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786
pub(super) fn tabular_spec(outfmt: &str) -> Option<TabularSpec> {
    let choice = crate::blastinput::app::parse_formatting_string(outfmt).ok()?;
    TabularSpec::of_app(
        crate::blastinput::app::ReportFormat::from_number(choice.number)?,
        &choice,
        crate::report::archive::options::RequestKind::Tblastn,
    )
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
    writer: &mut crate::cli::NcbiOstream<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    mut probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    ensure!(
        tab.queries.len() == hitlists.len() && tab.query_batch_skipped.len() == hitlists.len(),
        "TBLASTN tabular report query count mismatch"
    );
    let spec = tab.spec;
    // Which rows read their subject (`TabularSpec::subject_fetch`).
    let fetch = spec.subject_fetch(crate::blastinput::app::Program::Tblastn);
    let needs = spec.needs();
    let nident_alone = nident_alone(spec);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
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
        ShownSubjects::Fasta { parsed, .. } => TabularSubjects::Records {
            ids: tab.subject_ids,
            length: &record_length,
            parsed,
            records: tab.subject_records,
        },
        ShownSubjects::Database {
            reported,
            oids,
            seen,
            ..
        } => TabularSubjects::Database {
            ids: tab.subject_ids,
            ordinals: &reported.ordinals,
            oids: *oids,
            records: &reported.records,
            parsed: &reported.parsed,
            sources: &reported.sources,
            lengths: &lengths,
            taxonomy: &tab.inputs.taxonomy,
            loads: &reported.loads,
            seen: *seen,
        },
    };
    let names = tab.inputs.names.as_ref();
    // Outfmt 18, the organism report: the `Query=` block shows the query record's length, also
    // with `-query_loc` (as outfmt 0's).
    if spec.kind == TabularKind::TaxReport {
        let titles: Vec<_> = tab.queries.iter().map(FastaRecord::query_line).collect();
        let queries: Vec<crate::report::tax_report::TaxQuery<'_>> = tab
            .queries
            .iter()
            .zip(hitlists)
            .enumerate()
            .map(
                |(q_idx, (query, hitlist))| crate::report::tax_report::TaxQuery {
                    title: &titles[q_idx],
                    length: tab.ranges.queries.length(q_idx, query.seq().len()),
                    // A query of a batch that NCBI did not search has no Seq-align-set.
                    hsps: (!tab.query_batch_skipped[q_idx]).then(|| {
                        hitlist
                            .lists()
                            .iter()
                            .flat_map(|list| {
                                list.hsps.hsps.iter().zip(&list.payloads).map(
                                    move |(linked, payload)| crate::report::tax_report::TaxHsp {
                                        subject: list.oid as usize,
                                        evalue: linked.evalue,
                                        bit_score: payload.bit_score,
                                    },
                                )
                            })
                            .collect()
                    }),
                },
            )
            .collect();
        return crate::algorithm::common::tabular::write_tax_report(
            &subjects,
            tab.subject_records,
            false,
            names,
            PROGRAM,
            PROGRAM,
            &queries,
            writer,
            warnings,
            tab.inputs.is_html,
        );
    }
    // The refusals of `sall*` (`ParsedSubject::all_tabular_ids`) come before the first byte;
    // those of the other ids are `ReportedSubjects::tabular_subject_ids`'.
    if needs.all_deflines {
        for (q_idx, hitlist) in hitlists.iter().enumerate() {
            for rank in reported_ranks(std::slice::from_ref(hitlist)) {
                subjects.values_for(q_idx, rank, spec, PROGRAM)?;
            }
        }
    }
    let code = GeneticCode::try_from_id(tab.genetic_code).map_err(anyhow::Error::msg)?;
    // The positives of the rows: the formatter's `m_ScoringMatrix` (`GetAsciiProteinMatrix`).
    let table = TabularWriter::new(spec, tab.program_version, Some(tab.matrix));
    // NCBI makes a formatter for each query's report and for outfmt 7's epilog, and each warns
    // when the field list asks for a name of the taxonomy database that it cannot find.
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:139-158
    let taxdb_warning = spec.formatter_warning(tab.inputs.taxdb_found, "tblastn");
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
    // The bytes of one comment block or one row, written to the report stream as soon as they
    // are made (as NCBI's insertions reach its stream): nothing of a query's report is held.
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
        let query_length = tab.inputs.range_length.unwrap_or(searched_length);
        let tabular_query = TabularQuery::of_ids(query.believed(), query.shown_id(), query_length);
        let hit_count: usize = hitlist
            .lists()
            .iter()
            .map(|list| list.hsps.hsps.len())
            .sum();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180,204-208: a skipped search has a
        // null alignment set, so outfmt 7 has no count line (`write_tabular`).
        table.query_header(
            &mut buf,
            &query.query_line(),
            tab.shown.name().as_bytes(),
            (!tab.query_batch_skipped[q_idx]).then_some(hit_count),
        )?;
        writer.write_all(&buf)?;
        buf.clear();
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:816-824
        // `qcovs` reads the query range of every HSP of the query before its first row: a query
        // row that starts with a gap throws there, after outfmt 7's lines (`-subject`'s mapping
        // drops at most one end segment, `KappaResultHitList::put_in_report_order`, so any
        // remaining leading query gap still throws).
        if needs.qcovs {
            let unranged = hitlist
                .lists()
                .iter()
                .flat_map(|list| list.hsps.hsps.iter().zip(&list.payloads))
                .find_map(|(linked, payload)| {
                    reported_std_seg_row_range(
                        i64::from(linked.hsp.q_start),
                        i64::from(linked.hsp.q_end),
                        &payload.edit_script,
                        AlignRow::Query,
                        false,
                    )
                });
            if let Some(error) = unranged {
                writer.flush()?;
                return Err(NativeError::exception(&row_range_exception(error)).into());
            }
        }
        // The query's formatter (`TabularRows`).
        let mut rows = TabularRows::default();
        for list in hitlist.lists() {
            if list.hsps.hsps.is_empty() {
                continue;
            }
            let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
            ensure!(
                list.payloads.len() == list.hsps.hsps.len(),
                "TBLASTN report payload count mismatch"
            );
            let values = subjects.values_for(q_idx, oid, spec, PROGRAM)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:595-605
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
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:917,1061-1067
            // When NCBI's scope resolves the subject to a protein query Bioseq, `subject_is_na` is
            // false, so `sstrand` prints `N/A`.
            let seen_bioseq = subjects.seen_at(q_idx, oid);
            let mut subject = values.nucleotide_subject(names, qcovs, -1);
            if seen_bioseq.is_some() {
                subject.nucleotide = false;
            }
            // A subject whose lookups fail shows what the last found subject set
            // (`TabularRows::last_found`), made again from its index.
            let kept_values = match rows.last_found() {
                Some(kept) if table.fails(&subject) => {
                    Some(subjects.values_for(q_idx, kept, spec, PROGRAM)?)
                }
                _ => None,
            };
            let kept = kept_values
                .as_ref()
                .map(|values| values.nucleotide_subject(names, -1, -1));
            let subject_len = i32::try_from(tab.letters.length(oid))?;
            let whole_length = i32::try_from(record_length(oid))?;
            // The subject's `qcovhsp` rule: its rows only divide.
            let qcovhsp = QueryCoveragePerHsp::of_subject(
                tab.inputs.ungapped,
                list.hsps.hsps.len(),
                i32::try_from(searched_length)?,
                i32::try_from(query_length)?,
            );
            let mut translation = None;
            for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                if let Some(coding) = tab.query_other_coding.filter(|_| spec.reads_sequences()) {
                    return Err(crate::report::tabular::query_sequence_end(
                        writer, &buf, coding,
                    ));
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
                // A row that reads the subject (`TabularSpec::subject_fetch`, by the Seq-align's
                // identities for `nident` alone, `database::reads_subject`) loads it first where
                // the query's batch did not prefetch it: the texts of its taxonomy lookups follow
                // the rows before it (`QueryWarnings::subject_read`), before those of the row's
                // own lookups. A row that reads a subject whose load throws outside the subject
                // block ends the run (`TabularWriter::ends_run`).
                let engine_ident = usize::try_from(payload.num_ident).unwrap_or(0);
                if let Some(warnings) = warnings.as_deref_mut() {
                    if warnings.reads_with_texts(q_idx) && fetch.reads(engine_ident) {
                        warnings.subject_read(q_idx, oid as u32, &mut *writer)?;
                    }
                }
                if let Some(end) = table.ends_run(
                    &subject,
                    crate::blastinput::app::Program::Tblastn,
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
                let hsp = &linked.hsp;
                let frame = i32::from(hsp.frame);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1087-1105: the nucleotide
                // coordinates of a translated frame, descending on a minus frame (the default
                // row's, `write_tabular`).
                let (s_start, s_end) = super::stage_e_report::subject_endpoints(hsp, subject_len);
                let (q_start, q_end, s_start, s_end) = record_coordinates(
                    tab.ranges,
                    q_idx,
                    oid,
                    (hsp.q_start + 1, hsp.q_end, s_start, s_end),
                )?;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908
                // `num_ident` there is the Seq-align's `num_ident` score (the engine's identities,
                // `payload.num_ident`), which can differ from the report's identities. With
                // `nident` and no other field of the block, NCBI counts the identities only
                // when the engine's are positive and prints 0 otherwise. The shared writer
                // decides with `TabularHsp::num_ident` and counts again from the rows when
                // they are given: where the engine's identities are positive and the report's
                // are 0, the rows are read and `num_ident` is the engine's.
                // A row whose subject's load throws prints none of its letters.
                let rows_read = !subject.throws()
                    && (needs.rows
                        || (nident_alone && engine_ident > 0 && payload.report_num_ident == 0));
                if rows_read {
                    let translation = translation
                        .get_or_insert_with(|| TargetTranslation::of(tab.letters.nt(oid), &code));
                    // The rows without the lower case of the masks (`CAlnVec` of the query's
                    // Bioseq, whose letters the FASTA reader made upper case).
                    aligned_sequences_into(
                        QueryLetters::new(query_letters, tab.query_bioseq),
                        tab.letters.nt(oid),
                        &mut x_window,
                        &code,
                        translation,
                        hsp.frame,
                        hsp.q_start,
                        hsp.s_start,
                        hsp.s_end,
                        &payload.edit_script,
                        &mut qrow,
                        &mut srow,
                    )?;
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1072-1084
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1233-1244
                // The protein query has no strand (frame 0); the subject's frame is read in its
                // whole record, as the XML's `Hit_frame` (`translation_frame`).
                let s_frame = match seen_bioseq {
                    Some(bioseq) => {
                        crate::report::query_bioseq::ncbi_frame(frame > 0, s_start, bioseq.length())
                    }
                    None => translation_frame(
                        frame > 0,
                        s_start.min(s_end),
                        s_start.max(s_end),
                        whole_length,
                    ),
                };
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
                    ops: &payload.edit_script,
                    // The protein query's residues, in the HSP's context (with the letters that
                    // `-subject`'s mapper left out, which NCBI's BlastHSP still has) and in the
                    // Seq-align.
                    qcovhsp: qcovhsp.of_hsp(
                        hsp.q_end - hsp.q_start + payload.query_left_out as i32,
                        hsp.q_end - hsp.q_start,
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
                hsp_index += 1;
            }
        }
        // `~CBlastTabularInfo` of the query's report flushes the stream.
        crate::report::ostream::ReportWrite::ncbi_flush(writer)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
    // A search that stops at a query batch writes no epilog (`write_tabular`).
    if spec.kind == TabularKind::Commented && tab.ranges.epilog {
        warn(&mut *writer, warnings)?;
        table.epilog(writer, tab.queries.len())?;
        crate::report::ostream::ReportWrite::ncbi_flush(writer)?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
