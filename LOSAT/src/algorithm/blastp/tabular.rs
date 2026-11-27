//! BLASTP's (and psiblast's) tabular reports with a field list or a delimiter (outfmt 6 and 7
//! other than the default `std` fields with a tab) and outfmt 10 and 20: the adapter of the shared
//! writer (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query in order it writes the
//! query's warnings, outfmt 7's comment lines (psiblast's with its round) and one row per HSP of
//! the query's final hit list, in the list's order; outfmt 20's header line comes before the
//! first query of the run and outfmt 7's count of queries after the last.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:811-834
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
use super::blast_engine::query_rows::QueryRows;
use super::hsp::{BlastpHitList, BlastpHsp};
use crate::algorithm::common::tabular::FailedValues;
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::config::ScoringMatrix;
use crate::report::db_subject::ParsedSubject;
use crate::report::pairwise::PairwiseHit;
use crate::report::query_warnings::QueryWarnings;
use crate::report::tabular::{
    coverage_percent, union_length, DeflineIds, QueryCoveragePerHsp, SubjectLookup, TabularField,
    TabularHsp, TabularKind, TabularNames, TabularQuery, TabularRows, TabularSpec, TabularSubject,
    TabularTaxonomy, TabularWriter,
};

/// What the report needs besides the final hit lists (`Copy`: PHI-BLAST's report writes one view
/// of it for each query, `report_write::write_custom_tabular`).
#[derive(Clone, Copy)]
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
    /// The query records (the ids `-parse_deflines` read).
    pub query_records: &'a [FastaRecord],
    pub query_titles: &'a [Arc<[u8]>],
    pub searched_lengths: &'a [usize],
    /// The length of the `-query_loc` range as given (NCBI's `m_QueryRange`), `None` without it.
    pub range_length: Option<usize>,
    /// An ungapped search (`-ungapped`), whose Seq-aligns are Dense-diags (`qcovhsp`,
    /// `report::tabular::QueryCoveragePerHsp`).
    pub ungapped: bool,
    /// Outfmt 7's `# Database:` text.
    pub database: &'a [u8],
    pub subjects: TabularSubjects<'a>,
    /// By query: whether NCBI did not search its batch (outfmt 7's lines without a count).
    pub unsearched: &'a [bool],
    /// psiblast's count of outfmt 7's `# N hits found` line by query (the HSPs of its larger
    /// traceback hit list), `None` for blastp (the HSPs of the final hit list).
    pub hits_found: Option<&'a [usize]>,
    /// Whether outfmt 20's header line comes before the first query (the run's first report).
    pub prolog: bool,
    /// The number of queries formatted (outfmt 7's last line), `None` without the epilog.
    pub epilog: Option<usize>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath`).
    pub taxdb_found: bool,
    /// The names of the taxonomy database that the rows look up (`None`: every lookup fails).
    pub names: Option<&'a TabularNames>,
    /// The Seq-aligns' `num_ident` of the HSPs whose counts the reports' recount replaced
    /// (`recount_shown_identities`): NCBI's formatter opens its alignment block, and reads the
    /// subject for `nident`, on the Seq-align's count (`TabularSpec::opens_counts`).
    pub engine_counts: &'a super::seqalign::EngineCounts,
    /// The query's molecule is `other` (`PsiSearch::query_other_coding`): the first row that reads
    /// the query's letters ends the run (`report::tabular::query_sequence_end`).
    pub query_other_coding: Option<i32>,
}

/// The subjects as the report shows them, by subject index (a subject set) or rank (the
/// reported sequences of a database).
#[derive(Clone, Copy)]
pub(super) enum TabularSubjects<'a> {
    /// The records of `-subject`: the tabular id of each (`sseqid`, `sacc`, `saccver`; the
    /// default report's), the length of its record, and of the records whose ids `-parse_deflines`
    /// read their `sacc` and `saccver` and ids (both empty without the option).
    Records {
        records: &'a [FastaRecord],
        ids: &'a [crate::algorithm::common::tabular::SubjectId],
        accessions: &'a [Option<[Arc<[u8]>; 2]>],
        parsed: &'a [Option<Arc<ParsedSubject>>],
        length: &'a dyn Fn(usize) -> usize,
    },
    /// The reported sequences of a database: their tabular ids (`sseqid`, and `sacc` and
    /// `saccver` of a sequence with parsed Seq-ids; the default report's), ordinal ids, stored
    /// titles (the title of the first definition line of a sequence with parsed Seq-ids), parsed
    /// Seq-ids (empty when none has them), the sources from which NCBI makes the empty stored
    /// titles of the others (empty when none has one), lengths, and taxonomy (empty when no field
    /// shows it).
    Database {
        ids: &'a [crate::algorithm::common::tabular::SubjectId],
        accessions: &'a [Option<[Arc<[u8]>; 2]>],
        ordinals: &'a [u32],
        /// The OID of each ordinal id, which names a sequence without parsed Seq-ids.
        oids: crate::seqdb::OidMap<'a>,
        titles: &'a [Arc<[u8]>],
        parsed: &'a [Option<Arc<ParsedSubject>>],
        sources: &'a [crate::report::defline::SourceTaxname],
        length: &'a dyn Fn(usize) -> usize,
        taxonomy: &'a [TabularTaxonomy],
        /// What the reports do at each subject, by rank (`report::db_subject::HitLoads`): a
        /// subject whose load throws writes its rows as NCBI's `catch` leaves them
        /// (`report::tabular::SubjectLookup::Throws`).
        loads: &'a crate::report::db_subject::HitLoads,
        /// The queries whose Bioseqs NCBI's scope finds for reported sequences, by the call's
        /// query and the ordinal id (`report::query_bioseq`): such a row shows the query's
        /// (`TabularSubjects::values_for`).
        seen: crate::report::query_bioseq::SeenFrom<'a>,
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
    /// How the lookups of a subject fail (an empty id: `SetSubjectId`), `None` where they do not.
    failed: Option<FailedValues<'a>>,
}

impl<'a> SubjectValues<'a> {
    /// The values of the query Bioseq `bioseq` that NCBI's scope finds for a reported sequence
    /// (`common::tabular::SubjectValues::of_query_bioseq`).
    fn of_query_bioseq(
        bioseq: &'a crate::report::query_bioseq::QueryBioseq,
        all_deflines: bool,
    ) -> Self {
        let [id, accession, accession_version] = bioseq.tabular_ids();
        SubjectValues {
            id: id.as_slice(),
            accession: accession.as_slice(),
            accession_version: accession_version.as_slice(),
            gi: bioseq.gi(),
            all_ids: if all_deflines {
                bioseq.all_ids().to_vec()
            } else {
                Vec::new()
            },
            titles: None,
            length: bioseq.length(),
            taxonomy: Cow::Owned(bioseq.taxonomy()),
            failed: None,
        }
    }

    /// The values of the database sequence `subject` whose load throws with `exception`
    /// (`common::tabular::SubjectValues::thrown`): its rows print none of its own values, and
    /// NCBI reads none of its ids, titles or taxonomy.
    fn thrown(
        subject: &ParsedSubject,
        exception: &'a crate::ncbi_exception::NcbiException,
    ) -> Self {
        SubjectValues {
            id: &[],
            accession: &[],
            accession_version: &[],
            gi: 0,
            all_ids: Vec::new(),
            titles: None,
            length: 0,
            taxonomy: Cow::Owned(TabularTaxonomy::default()),
            failed: Some(FailedValues::thrown(subject, exception)),
        }
    }
}

impl SubjectValues<'_> {
    /// The subject of the shared writer, with the names of the taxonomy database (`None`: every
    /// lookup fails), `qcovs` of its HSPs (-1 for `N/A`) and its subject index `index`.
    fn subject<'s>(
        &'s self,
        names: Option<&'s TabularNames>,
        qcovs: i32,
        index: usize,
    ) -> TabularSubject<'s> {
        TabularSubject {
            id: self.id,
            accession: self.accession,
            accession_version: self.accession_version,
            gi: self.gi,
            all_ids: &self.all_ids,
            length: self.length,
            titles: self.titles.as_deref(),
            taxonomy: &self.taxonomy,
            names,
            nucleotide: false,
            qcovs,
            qcovus: -1,
            lookup: self.lookup(),
            index,
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

impl<'a> TabularSubjects<'a> {
    /// psiblast `-phi_pattern`: these subjects as the reports of the call's query `query`'s
    /// `results` pattern occurrences show them, each occurrence a result of its own
    /// (`report::query_bioseq::SeenFrom::of_occurrences`, whose value `occurrences` holds).
    #[cold]
    pub(crate) fn of_occurrences<'b>(
        &self,
        query: usize,
        results: usize,
        occurrences: &'b mut crate::report::query_bioseq::QueryBioseqs,
    ) -> TabularSubjects<'b>
    where
        'a: 'b,
    {
        match *self {
            TabularSubjects::Records {
                records,
                ids,
                accessions,
                parsed,
                length,
            } => TabularSubjects::Records {
                records,
                ids,
                accessions,
                parsed,
                length,
            },
            TabularSubjects::Database {
                ids,
                accessions,
                ordinals,
                oids,
                titles,
                parsed,
                sources,
                length,
                taxonomy,
                loads,
                seen,
            } => {
                *occurrences = seen.of_occurrences(query, results);
                TabularSubjects::Database {
                    ids,
                    accessions,
                    ordinals,
                    oids,
                    titles,
                    parsed,
                    sources,
                    length,
                    taxonomy,
                    loads,
                    seen: occurrences.numbered_from(0),
                }
            }
        }
    }

    /// The query Bioseq that NCBI's scope finds for subject `s_idx` in the report of the call's
    /// query `q_idx` (`report::query_bioseq::SeenFrom::of`; none for a subject set).
    pub(super) fn seen_at(
        &self,
        q_idx: usize,
        s_idx: usize,
    ) -> Option<&'a crate::report::query_bioseq::QueryBioseq> {
        match *self {
            TabularSubjects::Records { .. } => None,
            TabularSubjects::Database { ordinals, seen, .. } => {
                if !seen.any() {
                    return None;
                }
                seen.of(q_idx, *ordinals.get(s_idx)?)
            }
        }
    }

    /// `values` of subject `s_idx` in the report of the call's query `q_idx`: the values of the
    /// query Bioseq that NCBI's scope finds for it there (`seen_at`), else its own.
    fn values_for(
        &self,
        q_idx: usize,
        s_idx: usize,
        spec: &TabularSpec,
        app: BlastpApp,
    ) -> Result<SubjectValues<'a>> {
        match self.seen_at(q_idx, s_idx) {
            Some(bioseq) => Ok(SubjectValues::of_query_bioseq(
                bioseq,
                spec.needs().all_deflines,
            )),
            None => self.values(s_idx, spec, app),
        }
    }

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
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:560-565
            // A FASTA subject has no definition lines: its `sall*` are its own ids after
            // `s_ReplaceLocalId` (the title's first word, `FastaRecord::shown_id`, not the word of
            // `GetSeqIdList`), its title and taxonomy are empty (`N/A`).
            TabularSubjects::Records {
                records,
                ids,
                accessions,
                parsed,
                length,
            } => {
                // `None`: the lookup throws (`common::tabular::SubjectId`).
                let found = ids.get(s_idx).and_then(Option::as_deref);
                let id: &[u8] = found.unwrap_or_default();
                if let Some(subject) = parsed.get(s_idx).and_then(Option::as_deref) {
                    let (accession, accession_version): (&[u8], &[u8]) =
                        match accessions.get(s_idx).and_then(Option::as_ref) {
                            Some([accession, version]) => (accession, version),
                            None => (id, id),
                        };
                    let all = if needs.all_deflines {
                        subject.all_tabular_ids(true).map_err(|problem| {
                            anyhow::anyhow!(
                                "the subject record {} has {problem}, which is not supported by \
                                 LOSAT's {}",
                                String::from_utf8_lossy(&subject.alignment_fasta()),
                                app.name()
                            )
                        })?
                    } else {
                        Default::default()
                    };
                    return Ok(SubjectValues {
                        id,
                        accession,
                        accession_version,
                        gi: subject.tabular_gi(),
                        titles: None,
                        length: if shows_length { length(s_idx) } else { 0 },
                        taxonomy: Cow::Owned(TabularTaxonomy::default()),
                        failed: FailedValues::of_parsed(subject, found.is_none(), all.failed),
                        all_ids: all.ids,
                    });
                }
                // A record without a title: its `sseqid` is the first word of its generated title
                // (`unnamed`, `report::defline::tabular_subject_id`), its `sall*` its local id.
                let record = records.get(s_idx);
                let all_ids = match record {
                    // Only where the word is not the id, which the writer prints for none.
                    Some(record) if needs.all_deflines && record.shown_id() != id => {
                        DeflineIds::of_record(record.shown_id())
                    }
                    _ => Vec::new(),
                };
                // No id: `GetSeqIdList` throws (`FastaRecord::subject_tabular_id`).
                let failed = match record {
                    Some(record) if found.is_none() => {
                        let local = Cow::Borrowed(record.local_id.as_bytes());
                        Some(FailedValues {
                            align_id: [local.clone(), local],
                            subject_id: true,
                            all_ids: false,
                            thrown: None,
                            empty_handle: false,
                        })
                    }
                    _ => None,
                };
                Ok(SubjectValues {
                    id,
                    accession: id,
                    accession_version: id,
                    gi: 0,
                    all_ids,
                    titles: None,
                    length: if shows_length { length(s_idx) } else { 0 },
                    taxonomy: Cow::Owned(TabularTaxonomy::default()),
                    failed,
                })
            }
            TabularSubjects::Database {
                ids,
                accessions,
                ordinals,
                oids,
                titles,
                parsed,
                length,
                taxonomy,
                loads,
                ..
            } => {
                if let (Some(exception), Some(subject)) = (
                    loads.thrown(s_idx),
                    parsed.get(s_idx).and_then(Option::as_deref),
                ) {
                    return Ok(SubjectValues::thrown(subject, exception));
                }
                // `None`: the lookup throws (`common::tabular::SubjectId`).
                let found = ids.get(s_idx).and_then(Option::as_deref);
                let id: &[u8] = found.unwrap_or_default();
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
                        let all = if needs.all_deflines {
                            subject.all_tabular_ids(true).map_err(|problem| {
                                crate::algorithm::common::database::parsed_refusal(
                                    oids.oid(ordinals.get(s_idx).copied().unwrap_or_default()),
                                    &format!("has {problem}"),
                                    app.name(),
                                )
                            })?
                        } else {
                            Default::default()
                        };
                        let failed = FailedValues::of_parsed(subject, found.is_none(), all.failed);
                        Ok(SubjectValues {
                            id,
                            accession,
                            accession_version,
                            gi: subject.tabular_gi(),
                            all_ids: all.ids,
                            titles: Some(
                                subject
                                    .deflines()
                                    .iter()
                                    .map(|defline| defline.title.as_deref().unwrap_or_default())
                                    .collect(),
                            ),
                            length,
                            taxonomy,
                            failed,
                        })
                    }
                    // A sequence without parsed Seq-ids: its one definition line's title; no id,
                    // where `GetSeqIdList` throws for the general id
                    // (`TabularSubjectId::Unknown`), in `SetSubjectId` and `x_SetSubjectIds`.
                    None => Ok(SubjectValues {
                        id,
                        accession,
                        accession_version,
                        gi: 0,
                        all_ids: Vec::new(),
                        titles: Some(vec![titles.get(s_idx).map_or(&[][..], |title| title)]),
                        length,
                        taxonomy,
                        failed: found.is_none().then(|| {
                            let [full, accession] =
                                crate::algorithm::common::database::failed_general_lookup(
                                    oids.oid(ordinals.get(s_idx).copied().unwrap_or_default()),
                                );
                            FailedValues {
                                align_id: [Cow::Owned(full), Cow::Owned(accession)],
                                subject_id: true,
                                all_ids: true,
                                thrown: None,
                                empty_handle: false,
                            }
                        }),
                    }),
                }
            }
        }
    }
}

/// Whether `app`'s tabular report of `spec` reads the alignment rows (`TabularHsp::rows`): the
/// fields that show them (`TabularNeeds::rows`), and rpsblast's `nident`, `pident` and `mismatch`.
/// NCBI's formatter counts those on the aligned sequences of its scope, inside its alignment block
/// (`TabularSpec::opens_counts`), while the RPS engine's count, the Seq-align's `num_ident` that
/// opens the block for `nident`, is on the masked query (`core/blast_traceback.c:589-594`, the
/// core's subject is the user query as searched): rpsblast's rows of such a report are rendered
/// and handed to the HSPs whose block opens. The other programs count on the letters as read in
/// their engines, so the rows and the engine agree.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-980
/// `rows` is `TabularNeeds::rows` of `spec`.
pub(super) fn needs_rows(rows: bool, spec: &TabularSpec, app: BlastpApp) -> bool {
    rows || (app == BlastpApp::Rpsblast
        && [
            TabularField::NumIdentical,
            TabularField::PercentIdentical,
            TabularField::Mismatches,
        ]
        .into_iter()
        .any(|field| spec.has(field)))
}

/// The spec of a tabular report with a field list or a delimiter, or of outfmt 10 or 20, of a
/// `-outfmt` value read as NCBI reads it (`app::parse_formatting_string`), as `app`'s formatter
/// writes it (`TabularSpec::of_app`: psiblast's and rpsblast's without the custom delimiter);
/// `None` for the other formats and for the default outfmt 6 and 7. Every report of the blastp
/// engine reads its spec here, from the value of the command line or the library's as given.
pub(super) fn custom_tabular_spec(
    outfmt: impl AsRef<[u8]>,
    app: BlastpApp,
) -> Option<&'static TabularSpec> {
    let choice = crate::blastinput::app::parse_formatting_string(outfmt).ok()?;
    let format = crate::blastinput::app::ReportFormat::from_number(choice.number)?;
    TabularSpec::of_app(format, &choice, app.request_kind()).map(TabularSpec::intern)
}

/// Writes the warning of the formatter that NCBI makes for outfmt 20's header line, which it
/// writes before it reads the first query batch (`PrintProlog`, blastp_app.cpp:253): the
/// taxonomy database's warning when the field list checks it and `taxdb.bti` is not found
/// (`x_CheckTaxDB`, `TabularSpec::checks_taxdb`). The line itself comes with the first report
/// (`write_tabular_prolog`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
pub(super) fn write_tabular_prolog_warning(
    spec: &TabularSpec,
    taxdb_found: bool,
    program: &str,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if let Some(line) = spec
        .formatter_warning(taxdb_found, program)
        .filter(|_| spec.kind == TabularKind::CsvWithHeader)
    {
        diagnostics.write_all(line.as_bytes())?;
        diagnostics.flush()?;
    }
    Ok(())
}

/// Writes outfmt 20's line of the field names (`PrintProlog`; nothing for the other kinds).
pub(super) fn write_tabular_prolog(spec: &TabularSpec, writer: &mut dyn Write) -> Result<()> {
    TabularWriter::new(spec, "", None).prolog(writer)?;
    Ok(())
}

/// The reports of a run without a query batch: outfmt 20's line of the field names
/// (`PrintProlog`), outfmt 7's count of 0 queries with the warning of the epilog's formatter
/// (`PrintEpilog`), nothing for outfmt 6 and 10.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2238
pub(super) fn write_tabular_without_queries(
    spec: &TabularSpec,
    taxdb_found: bool,
    program: &str,
    writer: &mut dyn crate::report::ostream::ReportWrite,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    write_tabular_prolog(spec, writer)?;
    if spec.kind == TabularKind::Commented {
        if let Some(line) = spec.formatter_warning(taxdb_found, program) {
            writer.flush()?;
            diagnostics.write_all(line.as_bytes())?;
            diagnostics.flush()?;
        }
        TabularWriter::new(spec, "", None).epilog(writer, 0)?;
        // The epilog's `CBlastTabularInfo` flushes the stream in its destructor.
        writer.ncbi_flush()?;
    }
    Ok(())
}

/// psiblast's epilog of a tabular report after its last round (`blast_engine::psi_round`): outfmt
/// 7's count of the rounds formatted, after the warning of the epilog's formatter, as
/// `write_blastp_tabular` writes it with `tab.epilog`; nothing for the other kinds.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236 (cited at
// `write_blastp_tabular`)
#[cold]
#[inline(never)]
pub(super) fn write_tabular_round_epilog(
    spec: &TabularSpec,
    taxdb_found: bool,
    program: &str,
    writer: &mut dyn crate::report::ostream::ReportWrite,
    diagnostics: &mut dyn Write,
    rounds: usize,
) -> Result<()> {
    if spec.kind != TabularKind::Commented {
        return Ok(());
    }
    if let Some(line) = spec.formatter_warning(taxdb_found, program) {
        writer.flush()?;
        diagnostics.write_all(line.as_bytes())?;
        diagnostics.flush()?;
    }
    TabularWriter::new(spec, "", None).epilog(writer, rounds)?;
    writer.ncbi_flush()?;
    Ok(())
}

/// Writes the report of every query of `hit_lists` (indexed by query) to `writer`, each query's
/// warnings (and the taxonomy database's warning of its formatter) before its lines, outfmt 20's
/// header line first when `tab.prolog` and outfmt 7's last line after the last query when
/// `tab.epilog`. `probe` marks each row (the format observer). `rows` gives the rendered
/// alignments of each query's HSPs, in the order of its lists, when a field shows them
/// (`TabularNeeds::rows`; rpsblast's identities, `needs_rows`).
pub(super) fn write_blastp_tabular(
    tab: &BlastpTabular<'_>,
    hit_lists: &[Option<BlastpHitList>],
    writer: &mut dyn crate::report::ostream::ReportWrite,
    warnings: Option<&mut QueryWarnings<'_>>,
    probe: Option<&mut FormatProbe<'_>>,
    rows: Option<&mut QueryRows<'_>>,
) -> Result<()> {
    if tab.app == BlastpApp::Rpsblast {
        return write_rpsblast_tabular(tab, hit_lists, writer, warnings, probe, rows);
    }
    write_tabular_rows::<false>(tab, hit_lists, writer, warnings, probe, rows)
}

/// Writes the organism report (outfmt 18, `report::tax_report`) of every query of `hit_lists`
/// (indexed by query) to `writer`, each query's warnings before its report: the caption, the
/// `Query=` block of its title (`tab.query_titles`) with the length of its record
/// (`record_lengths`), and the table of its reported HSPs in the hit lists' order. The report is
/// made before the batch's first byte (`tax_report::report_batch`), so that LOSAT's refusals come
/// before it; `tab.names` carries the run's taxonomy database; `html`: `-html`.
pub(super) fn write_blastp_tax_report(
    tab: &BlastpTabular<'_>,
    record_lengths: &[usize],
    hit_lists: &[Option<BlastpHitList>],
    writer: &mut crate::cli::NcbiOstream<'_>,
    warnings: Option<&mut QueryWarnings<'_>>,
    html: bool,
) -> Result<()> {
    use crate::report::tax_report::{report_batch, write_report, TaxQuery};
    let queries: Vec<TaxQuery<'_>> = hit_lists
        .iter()
        .enumerate()
        .map(|(q_idx, hit_list)| TaxQuery {
            title: tab
                .query_titles
                .get(q_idx)
                .map_or(&b""[..], |title| title.as_ref()),
            length: record_lengths.get(q_idx).copied().unwrap_or_default(),
            // A query of a batch that NCBI did not search has no Seq-align-set.
            hsps: (!tab.unsearched.get(q_idx).copied().unwrap_or(false))
                .then(|| tax_hsps(hit_list.as_ref())),
        })
        .collect();
    let (parsed, plain, loads) = tax_subjects(tab);
    let seen = |q_idx: usize, s_idx: usize| tab.subjects.seen_at(q_idx, s_idx);
    let report = report_batch(
        &queries,
        parsed,
        plain,
        true,
        loads,
        tab.names,
        tab.app.name(),
        tab.app.ncbi_name(),
        html,
        Some(&seen),
    )?;
    write_report(writer, warnings, &queries, &report)?;
    Ok(())
}

/// The organism report's HSPs of the final hit list `hit_list` (`TaxHsp`), in its order.
#[inline(always)]
fn tax_hsps(hit_list: Option<&BlastpHitList>) -> Vec<crate::report::tax_report::TaxHsp> {
    hit_list
        .into_iter()
        .flat_map(reported)
        .flat_map(|hsp_list| &hsp_list.hsps)
        .map(|hsp| crate::report::tax_report::TaxHsp {
            subject: hsp.s_idx as usize,
            evalue: hsp.e_value,
            bit_score: hsp.bit_score,
        })
        .collect()
}

/// What the organism report reads of the subjects of `tab` (`tax_report::report_batch`): their
/// parsed Seq-ids, each subject's plain values (`TaxPlain`) and the database's loads.
#[inline(always)]
#[allow(clippy::type_complexity)]
fn tax_subjects<'t>(
    tab: &'t BlastpTabular<'t>,
) -> (
    &'t [Option<Arc<ParsedSubject>>],
    impl Fn(usize) -> crate::report::tax_report::TaxPlain<'t> + 't,
    Option<&'t crate::report::db_subject::HitLoads>,
) {
    use crate::report::tax_report::TaxPlain;
    let parsed = match &tab.subjects {
        TabularSubjects::Records { parsed, .. } | TabularSubjects::Database { parsed, .. } => {
            *parsed
        }
    };
    let plain = move |s_idx: usize| match &tab.subjects {
        TabularSubjects::Records { records, .. } => TaxPlain::Record {
            local_id: records[s_idx].local_id.as_bytes(),
            title: &records[s_idx].title,
        },
        TabularSubjects::Database {
            ordinals,
            oids,
            titles,
            sources,
            taxonomy,
            ..
        } => TaxPlain::Stored {
            oid: oids.oid(ordinals.get(s_idx).copied().unwrap_or_default()),
            title: titles.get(s_idx).map_or(&[][..], |title| title),
            source: crate::report::pairwise::stored_source(sources, s_idx as u32),
            taxid: taxonomy.get(s_idx).map_or(0, |taxonomy| taxonomy.taxid),
        },
    };
    let loads = match &tab.subjects {
        TabularSubjects::Records { .. } => None,
        TabularSubjects::Database { loads, .. } => Some(*loads),
    };
    (parsed, plain, loads)
}

/// One batch of PHI-BLAST's organism report (`write_phi_tax_report`): its queries (the call's),
/// its query with letters (`searched`), its messages, each result's query (whose title and record
/// length head its report) and final hit list (one for each occurrence), and whether NCBI
/// searched it.
pub(super) struct PhiTaxBatch {
    pub(super) queries: std::ops::Range<usize>,
    pub(super) searched: usize,
    pub(super) messages: Vec<u8>,
    pub(super) results: Vec<(usize, Option<BlastpHitList>)>,
    pub(super) unsearched: bool,
}

/// psiblast `-phi_pattern`: the organism report (outfmt 18) of PHI-BLAST's batches `batches`
/// (`PrintPhiResult`): for each, the warnings of its queries (`QueryWarnings::before_query`), its
/// messages, the caption once, then the report of each occurrence in turn
/// (`TaxReport::write_result`) under the title and length of its first query; a sequence that
/// the scope finds as the query shows it in each (`seen_at` of the batch's query with letters).
/// The reports are made before the first byte (`tax_report::report_batch`). The texts of the
/// subjects' reads are not placed by query here (the warnings' sink takes the lookups' texts).
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2058-2065
#[cold]
pub(super) fn write_phi_tax_report(
    tab: &BlastpTabular<'_>,
    record_lengths: &[usize],
    batches: &[PhiTaxBatch],
    writer: &mut crate::cli::NcbiOstream<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    html: bool,
) -> Result<()> {
    use crate::report::tax_report::{report_batch, TaxQuery};
    let mut queries = Vec::new();
    let mut query_of = Vec::new();
    for batch in batches {
        for (q_idx, list) in &batch.results {
            queries.push(TaxQuery {
                title: tab
                    .query_titles
                    .get(*q_idx)
                    .map_or(&b""[..], |title| title.as_ref()),
                length: record_lengths.get(*q_idx).copied().unwrap_or_default(),
                hsps: (!batch.unsearched).then(|| tax_hsps(list.as_ref())),
            });
            query_of.push(batch.searched);
        }
    }
    let (parsed, plain, loads) = tax_subjects(tab);
    let seen = |result: usize, s_idx: usize| tab.subjects.seen_at(query_of[result], s_idx);
    let report = report_batch(
        &queries,
        parsed,
        plain,
        true,
        loads,
        tab.names,
        tab.app.name(),
        tab.app.ncbi_name(),
        html,
        Some(&seen),
    )?;
    let mut result = 0;
    for batch in batches {
        if let Some(warnings) = warnings.as_deref_mut() {
            for q_idx in batch.queries.clone() {
                warnings.before_query(q_idx, &mut *writer)?;
            }
        }
        writer.write_all(&batch.messages)?;
        report.write_caption(writer)?;
        for _ in &batch.results {
            let query = &queries[result];
            match warnings.as_deref_mut() {
                Some(warnings) => {
                    let (_, sink) = warnings.reads_and_sink();
                    report.write_result(result, writer, sink, None, query.title, query.length)?;
                }
                None => report.write_result(
                    result,
                    writer,
                    &mut std::io::sink(),
                    None,
                    query.title,
                    query.length,
                )?,
            }
            result += 1;
        }
    }
    if let Some(warnings) = warnings {
        warnings.after_last_query(batches.last().map_or(0, |batch| batch.queries.end), writer)?;
    }
    Ok(())
}

/// rpsblast's report (`write_tabular_rows` with the gate of its rows on NCBI's alignment block,
/// `needs_rows`). Cold and out of line (no NCBI counterpart: code placement only), so that
/// blastp's and psiblast's writer keeps the default path's code.
#[cold]
#[inline(never)]
fn write_rpsblast_tabular(
    tab: &BlastpTabular<'_>,
    hit_lists: &[Option<BlastpHitList>],
    writer: &mut dyn crate::report::ostream::ReportWrite,
    warnings: Option<&mut QueryWarnings<'_>>,
    probe: Option<&mut FormatProbe<'_>>,
    rows: Option<&mut QueryRows<'_>>,
) -> Result<()> {
    write_tabular_rows::<true>(tab, hit_lists, writer, warnings, probe, rows)
}

/// `write_blastp_tabular`; `RPS`: rpsblast's report, whose rows are read for its identities and
/// handed only to the HSPs whose alignment block opens (`needs_rows`).
fn write_tabular_rows<const RPS: bool>(
    tab: &BlastpTabular<'_>,
    hit_lists: &[Option<BlastpHitList>],
    writer: &mut dyn crate::report::ostream::ReportWrite,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut rows: Option<&mut QueryRows<'_>>,
) -> Result<()> {
    let spec = tab.spec;
    let needs = spec.needs();
    // Which rows read their subject (`TabularSpec::subject_fetch`).
    let fetch = spec.subject_fetch(crate::blastinput::app::Program::Blastp);
    let rows_read = if RPS {
        needs_rows(needs.rows, spec, tab.app)
    } else {
        needs.rows
    };
    // The refusals of `sall*` (`ParsedSubject::all_tabular_ids`) come before the first byte;
    // those of the other ids are made with them (`tabular_subject_ids`).
    if needs.all_deflines {
        for (q_idx, hit_list) in hit_lists.iter().enumerate() {
            for hsp_list in hit_list.iter().flat_map(reported) {
                if let Some(first) = hsp_list.hsps.first() {
                    tab.subjects
                        .values_for(q_idx, first.s_idx as usize, spec, tab.app)?;
                }
            }
        }
    }
    // rpsblast's options carry no matrix name (the `.aux` matrix goes to the search's options
    // only), so NCBI's formatter has no matrix and counts the identities as the positives.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:238
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3552-3556
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1017-1020
    let table = TabularWriter::new(
        spec,
        tab.program_version,
        (tab.app != BlastpApp::Rpsblast).then_some(tab.matrix),
    );
    // NCBI makes a formatter for each query's report and for outfmt 7's epilog, and each warns
    // when the field list asks for a name of the taxonomy database that it cannot find.
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:139-158
    let taxdb_warning = spec.formatter_warning(tab.taxdb_found, tab.app.ncbi_name());
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
    // The bytes of one comment block or one row, written to the report stream as soon as they
    // are made (as NCBI's insertions reach its stream): nothing of a query's report is held.
    let mut buf: Vec<u8> = Vec::new();
    let mut query_ranges: Vec<(u64, u64)> = Vec::new();
    // The index of the next HSP in the final hit lists (the probe's HSP number).
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
        let query_length = tab.range_length.unwrap_or(searched_length);
        let query = TabularQuery::of_ids(
            tab.query_records.get(q_idx).and_then(FastaRecord::believed),
            tab.query_ids.get(q_idx).map_or(&b""[..], |id| id.as_ref()),
            query_length,
        );
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1277-1278
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
            tab.database,
            num_hits,
        )?;
        writer.write_all(&buf)?;
        buf.clear();
        let Some(hit_list) = hit_list else {
            // `~CBlastTabularInfo` of the query's report flushes the stream.
            writer.ncbi_flush()?;
            continue;
        };
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-986
        // The rows of the query's HSPs when a field shows them, and the next one's index.
        let query_rows: Option<&[PairwiseHit]> = match rows.as_deref_mut() {
            Some(rows) if rows_read => Some(rows.query(q_idx)?),
            _ => None,
        };
        let mut row_index = 0usize;
        // The query's (or round's) formatter (`TabularRows`).
        let mut table_rows = TabularRows::default();
        for (list_index, hsp_list) in reported(hit_list).enumerate() {
            let Some(first) = hsp_list.hsps.first() else {
                continue;
            };
            let s_idx = first.s_idx as usize;
            let values = tab.subjects.values_for(q_idx, s_idx, spec, tab.app)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:595-605
            // The HSPs of one subject are one group of `InsertSubjectScores` (a hit list holds a
            // subject once); `qcovus` is blastn's only (blast_format.cpp:818-820).
            let qcovs = if needs.qcovs {
                query_ranges.clear();
                // An HSP of length 0 has no letter to cover (its range is empty).
                query_ranges.extend(
                    hsp_list
                        .hsps
                        .iter()
                        .filter(|hsp| hsp.length > 0)
                        .map(|hsp| (hsp.q_start as u64 - 1, hsp.q_end as u64 - 1)),
                );
                coverage_percent(union_length(&mut query_ranges), query_length)
            } else {
                -1
            };
            let subject = values.subject(tab.names, qcovs, s_idx);
            // A subject whose lookups fail shows what the last found subject set
            // (`TabularRows::last_found`), made again from its index.
            let kept_values = match table_rows.last_found() {
                Some(kept) if table.fails(&subject) => {
                    Some(tab.subjects.values_for(q_idx, kept, spec, tab.app)?)
                }
                _ => None,
            };
            let kept = kept_values
                .as_ref()
                .zip(table_rows.last_found())
                .map(|(values, kept)| values.subject(tab.names, -1, kept));
            // The subject's `qcovhsp` rule: its rows only divide.
            let qcovhsp = QueryCoveragePerHsp::of_subject(
                tab.ungapped,
                hsp_list.hsps.len(),
                searched_length as i32,
                query_length as i32,
            );
            for (index_in_list, hsp) in hsp_list.hsps.iter().enumerate() {
                if let Some(coding) = tab.query_other_coding.filter(|_| spec.reads_sequences()) {
                    return Err(crate::report::tabular::query_sequence_end(
                        writer, &buf, coding,
                    ));
                }
                // The Seq-align's `num_ident` (the engine's count where the reports' recount
                // replaced it, `BlastpTabular::engine_counts`).
                let seqalign_num_ident = tab.engine_counts.num_ident(
                    (q_idx as u32, list_index as u32, index_in_list as u32),
                    hsp.num_ident,
                );
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
                // A row that reads the subject loads it first where the query's batch did not
                // prefetch it: the texts of its taxonomy lookups follow the rows before it
                // (`QueryWarnings::subject_read`), before those of the row's own lookups.
                if let Some(warnings) = warnings.as_deref_mut() {
                    if warnings.reads_with_texts(q_idx) && fetch.reads(seqalign_num_ident) {
                        warnings.subject_read(q_idx, hsp.s_idx, &mut *writer)?;
                    }
                }
                // A row that reads a subject whose load throws outside the subject block ends
                // the run (`TabularWriter::ends_run`).
                if let Some(exception) = table.ends_run(
                    &subject,
                    crate::blastinput::app::Program::Blastp,
                    seqalign_num_ident,
                ) {
                    return Err(crate::report::tabular::run_end(
                        writer,
                        &buf,
                        exception,
                        || {
                            if let Some(warnings) = warnings.as_deref_mut() {
                                warnings.end_run();
                            }
                        },
                    ));
                }
                // The texts of the row's lookups of the taxonomy database (`x_SetTaxInfoAll`,
                // `x_SetTaxInfo`) follow the rows before it.
                if let (Some(names), Some(warnings), true) = (
                    tab.names,
                    warnings.as_deref_mut(),
                    table.sets_subject(&subject),
                ) {
                    let text = names.row_messages(spec, &values.taxonomy);
                    if !text.is_empty() {
                        writer.flush()?;
                        warnings.sink.write_all(&text)?;
                        warnings.sink.flush()?;
                    }
                }
                let rows = match query_rows {
                    // rpsblast's identities: the rows inside NCBI's alignment block only (`nident`
                    // alone opens it on the engine's count, `needs_rows`).
                    Some(_) if RPS && !spec.opens_counts(&seqalign_num_ident) => None,
                    Some(rows) => Some(hsp_rows(rows, row_index, hsp)?),
                    None => None,
                };
                row_index += 1;
                let mut row = tabular_hsp(hsp, qcovhsp, rows);
                // Outside NCBI's alignment block `nident` is the 0 of `x_ResetFields`: with
                // `nident` the block stays closed only on a Seq-align count of 0, which is that 0.
                if !spec.opens_counts(&seqalign_num_ident) {
                    row.num_ident = seqalign_num_ident;
                }
                // One printed row is one HSP; the probe marks it without changing it.
                if let Some(probe) = probe.as_deref_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                table.row(
                    &mut buf,
                    &mut table_rows,
                    &query,
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
        writer.ncbi_flush()?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
    // The warnings after the last query's report (`QueryWarnings::after_last_query`).
    if let Some(warnings) = warnings.as_deref_mut() {
        warnings.after_last_query(hit_lists.len(), &mut *writer)?;
    }
    if let (TabularKind::Commented, Some(num_queries)) = (spec.kind, tab.epilog) {
        warn(&mut *writer, warnings)?;
        table.epilog(writer, num_queries)?;
        writer.ncbi_flush()?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:271-327
/// The HSP lists of a final hit list that the reports show (its first `hsplist_count`), as the
/// default report walks them (`write_blastp_tabular_hit_lists`).
fn reported(hit_list: &BlastpHitList) -> impl Iterator<Item = &super::hsp::BlastpHspList> {
    hit_list.hsplist_array.iter().take(hit_list.hsplist_count)
}

/// The rows of the HSP `hsp` (the `index`-th of its query's final list, whose rows are `rows`):
/// its rendered alignment, whose end gaps the table reads from the HSP's script
/// (`TabularHsp::ops`).
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
/// its `qcovhsp` by its subject's rule `coverage`, `rows` when a field shows them.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1071-1094
/// The gaps are the alignment's columns that are neither identities nor mismatches
/// (`BlastpHsp::gap_letters`); the gap openings and the identities are the default row's. A
/// protein HSP runs forward on both sequences: the frames are 1.
fn tabular_hsp<'r>(
    hsp: &'r BlastpHsp,
    coverage: QueryCoveragePerHsp,
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
        gaps: hsp.gap_letters(),
        num_ident: hsp.num_ident,
        q_frame: 1,
        s_frame: 1,
        rows,
        ops: hsp.gap_info.as_deref().unwrap_or_default(),
        // The HSP's query letters, in its context and in its Seq-align.
        qcovhsp: {
            let span = (hsp.q_end + 1 - hsp.q_start) as i32;
            coverage.of_hsp(span, span)
        },
    }
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
