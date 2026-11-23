//! The subjects of the tabular reports with a field list (outfmt 6 and 7 other than the default
//! `std` fields with a tab, outfmt 10 and 20; `report::tabular`), as the programs' adapters hand
//! them to the shared writer: the values of a FASTA subject or a database sequence that its rows
//! show, the taxonomy of the reported database sequences, and outfmt 20's line of the field names
//! with its formatter's warning. For the programs whose subjects are nucleotide sequences;
//! `program` names the program in LOSAT's refusals. BLASTN (line CLI C2, PR 1; its outfmt 20 line
//! and warning are its own, DECISIONS 6), TBLASTN and TBLASTX (PR 2) use it.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-851

use std::borrow::Cow;
use std::io::Write;
use std::sync::Arc;

use anyhow::Result;

use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::db_subject::ParsedSubject;
use crate::report::tabular::{
    taxdb_warning_line, DeflineIds, TabularField, TabularKind, TabularNames, TabularSpec,
    TabularSubject, TabularTaxonomy, TabularWriter,
};
use crate::seqdb::SeqDb;

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908
/// The fields that enter `SetFields`' block of the coordinates, the rows and the counts
/// (`SetCounts`), besides `nident` with a positive engine count (`num_ident` of the Seq-align).
pub(crate) const COUNTS_BLOCK: [TabularField; 16] = [
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

/// The tabular id of one subject that the default fields print (`saccver`; for a subject without
/// parsed Seq-ids also its `sseqid` and `sacc`), as the programs' lists of them hold it by subject
/// index or rank: `None` where NCBI's `GetSeqIdList` makes none and `SetSubjectId` throws
/// (`report::tabular::SubjectLookup::Failed`), whose rows NCBI writes with its catch's values or
/// stops at (`report::tabular::check_failed_lookups`).
pub(crate) type SubjectId = Option<Arc<[u8]>>;

/// The bytes of `ids[s_idx]` (`SubjectId`) that a row prints: empty where the lookup fails, whose
/// default fields' rows are refused before they are written (`check_failed_lookups`), or where
/// `ids` has no such subject.
#[inline]
pub(crate) fn printed_id(ids: &[SubjectId], s_idx: usize) -> &[u8] {
    ids.get(s_idx)
        .and_then(Option::as_deref)
        .unwrap_or_default()
}

/// The subjects as the report shows them, by subject index (a subject set) or rank (the
/// reported sequences of a database).
pub(crate) enum TabularSubjects<'a> {
    /// The records of `-subject`: the tabular id of each (`sseqid`, `sacc`, `saccver`; the
    /// default report's), the length of its record, and the ones whose ids `-parse_deflines`
    /// read (empty without the option).
    Records {
        ids: &'a [SubjectId],
        length: &'a dyn Fn(usize) -> usize,
        parsed: &'a [Option<Arc<ParsedSubject>>],
        /// The records, by subject index: the Seq-align id (`FastaRecord::local_id`) of one whose
        /// lookup fails (`None` in `ids`) and the first word of the title that `sall*` print.
        records: &'a [FastaRecord],
    },
    /// The reported sequences of a database: their tabular ids (`saccver`,
    /// `common::database::tabular_subject_ids`), ordinal ids, records (with
    /// the stored title of a sequence without parsed Seq-ids), parsed Seq-ids (empty when none
    /// has them), the names of the sources from which NCBI makes empty stored titles (empty when
    /// none has one), lengths, and taxonomy (empty when no field shows it).
    Database {
        ids: &'a [SubjectId],
        ordinals: &'a [u32],
        /// The OID of each ordinal id, which names a sequence without parsed Seq-ids.
        oids: crate::seqdb::OidMap<'a>,
        records: &'a [FastaRecord],
        parsed: &'a [Option<Arc<ParsedSubject>>],
        sources: &'a [crate::report::defline::SourceTaxname],
        lengths: &'a [usize],
        taxonomy: &'a [TabularTaxonomy],
        /// What the reports do at each subject (`report::db_subject::HitLoads`, by the subjects'
        /// key here): a subject whose load throws writes its rows as NCBI's `catch` leaves them
        /// (`report::tabular::SubjectLookup::Throws`).
        loads: &'a crate::report::db_subject::HitLoads,
    },
}

/// The values of one subject that its rows show (`TabularSubject` without the coverages).
pub(crate) struct SubjectValues<'a> {
    id: Cow<'a, [u8]>,
    accession: Cow<'a, [u8]>,
    accession_version: Cow<'a, [u8]>,
    gi: i64,
    all_ids: Vec<DeflineIds>,
    titles: Option<Vec<&'a [u8]>>,
    length: usize,
    taxonomy: Cow<'a, TabularTaxonomy>,
    /// How the lookups of a subject fail (`FailedLookup`), `None` where they do not.
    failed: Option<FailedValues<'a>>,
    /// The subject index (`TabularSubject::index`).
    index: usize,
}

/// The values of `FailedLookup` that a subject's rows show where its lookups fail.
pub(crate) struct FailedValues<'a> {
    /// The Seq-align id as `sallseqid` and `sallacc` print it (`FailedLookup::full`,
    /// `accession`).
    pub(crate) align_id: [Cow<'a, [u8]>; 2],
    /// Whether `SetSubjectId` throws (`FailedLookup::subject_id`).
    pub(crate) subject_id: bool,
    /// Whether `x_SetSubjectIds` throws, after the lines that `SubjectValues::all_ids` then holds
    /// (`FailedLookup::all_ids`).
    pub(crate) all_ids: bool,
    /// Where the load of the subject throws before both (`report::tabular::ThrownLookup`): the
    /// Seq-align id's GI and the exception.
    pub(crate) thrown: Option<(i64, &'a crate::ncbi_exception::NcbiException)>,
    /// Whether the scope's Bioseq handle of the subject is empty (`FailedLookup::empty_handle`).
    pub(crate) empty_handle: bool,
}

impl<'a> FailedValues<'a> {
    /// The failure of a `-subject` record or a database sequence with parsed Seq-ids whose
    /// `SetSubjectId` throws where `subject_id` (`ParsedSubject::tabular_ids` gave none) and whose
    /// `x_SetSubjectIds` throws where `all_ids` (`AllTabularIds::failed`), `None` where neither.
    /// A database sequence whose Bioseq the scope does not get (`ParsedSubject::scope_finds_none`)
    /// fails on its empty handle instead (`FailedLookup::empty_handle`): `SetSubjectId` throws,
    /// and `GetBioseqLength` before `x_SetSubjectIds`.
    pub(crate) fn of_parsed(
        subject: &ParsedSubject,
        subject_id: bool,
        all_ids: bool,
    ) -> Option<Self> {
        let empty_handle = subject.scope_finds_none();
        (subject_id || all_ids || empty_handle).then(|| {
            let [full, accession] = subject.alignment_tabular_ids();
            FailedValues {
                align_id: [Cow::Owned(full), Cow::Owned(accession)],
                subject_id: subject_id || empty_handle,
                all_ids: all_ids && !empty_handle,
                thrown: None,
                empty_handle,
            }
        })
    }

    /// The failure of a database sequence with parsed Seq-ids whose load throws with `exception`
    /// where the query's batch has its prefetch off (`report::db_subject::HitLoads::thrown`).
    pub(crate) fn thrown(
        subject: &ParsedSubject,
        exception: &'a crate::ncbi_exception::NcbiException,
    ) -> Self {
        let [full, accession] = subject.alignment_tabular_ids();
        FailedValues {
            align_id: [Cow::Owned(full), Cow::Owned(accession)],
            subject_id: false,
            all_ids: false,
            thrown: Some((subject.alignment_gi(), exception)),
            empty_handle: false,
        }
    }

    /// The lookup of a subject that fails as these values say, with the definition lines' ids
    /// `all_ids` that `x_SetSubjectIds` pushed before it threw (`TabularSubject::lookup`).
    pub(crate) fn lookup<'s>(
        &'s self,
        all_ids: &'s [crate::report::tabular::DeflineIds],
    ) -> crate::report::tabular::SubjectLookup<'s> {
        match self.thrown {
            Some((gi, exception)) => crate::report::tabular::SubjectLookup::Throws(
                crate::report::tabular::ThrownLookup {
                    full: &self.align_id[0],
                    gi,
                    accession: &self.align_id[1],
                    exception,
                },
            ),
            None => crate::report::tabular::SubjectLookup::Failed(
                crate::report::tabular::FailedLookup {
                    full: &self.align_id[0],
                    accession: &self.align_id[1],
                    subject_id: self.subject_id,
                    all_ids: self.all_ids.then_some(all_ids),
                    empty_handle: self.empty_handle,
                },
            ),
        }
    }
}

impl<'a> SubjectValues<'a> {
    /// The values of the database sequence `subject` (subject index `index`) whose load throws
    /// with `exception` (`FailedValues::thrown`): its rows print none of its own values (NCBI's
    /// `catch` leaves the subject block's members as the rows before set them), and NCBI reads
    /// none of its ids, titles or taxonomy, which LOSAT does not check.
    pub(crate) fn thrown(
        subject: &ParsedSubject,
        exception: &'a crate::ncbi_exception::NcbiException,
        index: usize,
    ) -> Self {
        SubjectValues {
            id: Cow::Borrowed(&[]),
            accession: Cow::Borrowed(&[]),
            accession_version: Cow::Borrowed(&[]),
            gi: 0,
            all_ids: Vec::new(),
            titles: None,
            length: 0,
            taxonomy: Cow::Owned(TabularTaxonomy::default()),
            failed: Some(FailedValues::thrown(subject, exception)),
            index,
        }
    }
}

impl SubjectValues<'_> {
    /// The subject of the shared writer: a nucleotide sequence (`sstrand` is `plus` or `minus`)
    /// with the names of the taxonomy database that its rows look up (`None`: every lookup
    /// fails) and the coverages `qcovs` and `qcovus` of its HSPs (-1 for `N/A`).
    pub(crate) fn nucleotide_subject<'s>(
        &'s self,
        names: Option<&'s TabularNames>,
        qcovs: i32,
        qcovus: i32,
    ) -> TabularSubject<'s> {
        TabularSubject {
            id: &self.id,
            accession: &self.accession,
            accession_version: &self.accession_version,
            gi: self.gi,
            all_ids: &self.all_ids,
            length: self.length,
            titles: self.titles.as_deref(),
            taxonomy: &self.taxonomy,
            names,
            nucleotide: true,
            qcovs,
            qcovus,
            lookup: self.lookup(),
            index: self.index,
        }
    }

    /// The lookup of the subject (`TabularSubject::lookup`).
    pub(crate) fn lookup(&self) -> crate::report::tabular::SubjectLookup<'_> {
        match &self.failed {
            Some(failed) => failed.lookup(&self.all_ids),
            None => crate::report::tabular::SubjectLookup::Found,
        }
    }
}

/// Which lookups of subject `s_idx` fail (`report::tabular::check_failed_lookups`) for a report of
/// `spec` (`None`: the default fields): its id among `ids` (`SubjectId`: `None` where
/// `SetSubjectId` throws), and for a subject with parsed Seq-ids (`parsed`) and a report of
/// `sall*`, also `x_SetSubjectIds` (`ParsedSubject::lookup_failure`). `database`: the subjects are
/// a database's (the one definition line of a sequence without parsed Seq-ids throws in
/// `x_SetSubjectIds` too, `LookupFailure::of_plain`); `protein` is their kind. A subject whose
/// load throws (`loads`, by subject index: `report::db_subject::HitLoad::Throws`) fails there,
/// before every other lookup (`LookupFailure::thrown`).
pub(crate) fn subject_lookup_failure(
    spec: Option<&TabularSpec>,
    ids: &[SubjectId],
    parsed: &[Option<Arc<ParsedSubject>>],
    loads: &crate::report::db_subject::HitLoads,
    database: bool,
    protein: bool,
    s_idx: usize,
) -> Option<crate::report::tabular::LookupFailure> {
    if loads.thrown(s_idx).is_some() {
        return Some(crate::report::tabular::LookupFailure::thrown());
    }
    let failed_id = matches!(ids.get(s_idx), Some(None));
    match parsed.get(s_idx).and_then(Option::as_deref) {
        Some(subject) if subject.scope_finds_none() => {
            Some(crate::report::tabular::LookupFailure::empty_handle())
        }
        Some(subject) if spec.is_some_and(reports_all_ids) => {
            subject.lookup_failure(protein).ok().flatten()
        }
        Some(_) => failed_id.then_some(crate::report::tabular::LookupFailure {
            subject_id: true,
            all_ids: false,
            bioseq: false,
            length: false,
        }),
        None => failed_id.then(|| crate::report::tabular::LookupFailure::of_plain(database)),
    }
}

/// Whether the report shows `sallseqid`, `sallgi` or `sallacc` (`x_SetSubjectIds` runs).
pub(crate) fn reports_all_ids(spec: &TabularSpec) -> bool {
    spec.has_any(&[
        TabularField::SubjectAllSeqIds,
        TabularField::SubjectAllGis,
        TabularField::SubjectAllAccessions,
    ])
}

/// Whether the lookups of some subject can fail, so that a tabular report may need the check of
/// `checks_failed_lookups`: a failed id among `ids` (`SetSubjectId` throws), or subjects with
/// parsed Seq-ids (`parsed` is empty when none has them).
#[inline]
pub(crate) fn lookups_may_fail(ids: &[SubjectId], parsed: &[Option<Arc<ParsedSubject>>]) -> bool {
    !parsed.is_empty() || ids.iter().any(Option::is_none)
}

/// Whether a tabular report of `spec` is checked for rows whose lookups fail
/// (`report::tabular::check_failed_lookups`) where the subjects' `ids` (`None` for a failed
/// `SetSubjectId`) and `parsed` are as given: some id failed, a subject with parsed Seq-ids
/// may fail in `x_SetSubjectIds` (`subject_lookup_failure`) or has no Bioseq
/// (`ParsedSubject::scope_finds_none`), or the load of some subject throws (`loads`).
pub(crate) fn checks_failed_lookups(
    spec: Option<&TabularSpec>,
    ids: &[SubjectId],
    parsed: &[Option<Arc<ParsedSubject>>],
    loads: &crate::report::db_subject::HitLoads,
) -> bool {
    (ids.iter().any(Option::is_none)
        || loads.any_throws()
        || (!parsed.is_empty()
            && spec.is_some_and(reports_all_ids)
            && parsed.iter().any(Option::is_some))
        || parsed
            .iter()
            .flatten()
            .any(|subject| subject.scope_finds_none()))
        && crate::report::tabular::may_stop_at_failed_lookup(spec)
}

impl<'a> TabularSubjects<'a> {
    /// The parsed subjects by subject index (empty when none has parsed Seq-ids).
    pub(crate) fn parsed(&self) -> &'a [Option<Arc<ParsedSubject>>] {
        match *self {
            TabularSubjects::Records { parsed, .. } | TabularSubjects::Database { parsed, .. } => {
                parsed
            }
        }
    }

    /// What the reports load of each subject (`report::db_subject::HitLoads`; `None` for a
    /// subject set, whose subjects are all loaded).
    pub(crate) fn loads(&self) -> Option<&'a crate::report::db_subject::HitLoads> {
        match *self {
            TabularSubjects::Records { .. } => None,
            TabularSubjects::Database { loads, .. } => Some(loads),
        }
    }

    /// Subject `s_idx` without its parsed Seq-ids as the organism report (outfmt 18,
    /// `report::tax_report`) reads it; `records` are the records of a subject set (a database's
    /// are its own).
    pub(crate) fn tax_plain(
        &self,
        s_idx: usize,
        records: &'a [FastaRecord],
    ) -> crate::report::tax_report::TaxPlain<'a> {
        use crate::report::tax_report::TaxPlain;
        match *self {
            TabularSubjects::Records { .. } => TaxPlain::Record {
                local_id: records[s_idx].local_id.as_bytes(),
                title: &records[s_idx].title,
            },
            TabularSubjects::Database {
                ordinals,
                oids,
                records,
                sources,
                taxonomy,
                ..
            } => TaxPlain::Stored {
                oid: oids.oid(ordinals.get(s_idx).copied().unwrap_or_default()),
                title: records.get(s_idx).map_or(&[][..], |record| &record.title),
                source: sources
                    .get(s_idx)
                    .unwrap_or(&crate::report::defline::UNNAMED_SOURCE),
                taxid: taxonomy.get(s_idx).map_or(0, |taxonomy| taxonomy.taxid),
            },
        }
    }

    /// The values of subject `s_idx` that the fields of `spec` show; `program` names the program
    /// in LOSAT's refusal of a sequence with parsed Seq-ids whose ids it does not write as NCBI.
    pub(crate) fn values(
        &self,
        s_idx: usize,
        spec: &TabularSpec,
        program: &str,
    ) -> Result<SubjectValues<'a>> {
        let needs = spec.needs();
        match *self {
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:560-565
            // A FASTA subject has no definition lines: its `sall*` are its own ids after
            // `s_ReplaceLocalId` (the title's first word, `FastaRecord::shown_id`, not the word of
            // `GetSeqIdList`), its title and taxonomy are empty (`N/A`).
            TabularSubjects::Records {
                ids,
                length,
                parsed,
                records,
            } => {
                let length = if spec.has(TabularField::SubjectLength) {
                    length(s_idx)
                } else {
                    0
                };
                if let Some(subject) = parsed.get(s_idx).and_then(Option::as_deref) {
                    let refused = |problem: &str| {
                        anyhow::anyhow!(
                            "the subject record {} has {problem}, which is not supported by \
                             LOSAT's {program}",
                            String::from_utf8_lossy(&subject.alignment_fasta())
                        )
                    };
                    let ids = subject.tabular_ids(false).map_err(refused)?;
                    let all = if needs.all_deflines {
                        subject.all_tabular_ids(false).map_err(refused)?
                    } else {
                        Default::default()
                    };
                    let failed = FailedValues::of_parsed(subject, ids.is_none(), all.failed);
                    let [id, accession, accession_version] = ids.unwrap_or_default();
                    return Ok(SubjectValues {
                        id: Cow::Owned(id),
                        accession: Cow::Owned(accession),
                        accession_version: Cow::Owned(accession_version),
                        gi: subject.tabular_gi(),
                        all_ids: all.ids,
                        titles: None,
                        length,
                        taxonomy: Cow::Owned(TabularTaxonomy::default()),
                        failed,
                        index: s_idx,
                    });
                }
                let found = ids.get(s_idx).map(Option::as_deref);
                let id: &[u8] = found.map_or(&b"unknown"[..], Option::unwrap_or_default);
                let record = records.get(s_idx);
                let all_ids = match record {
                    // Only where the word is not the id, which the writer prints for none.
                    Some(record) if needs.all_deflines && record.shown_id() != id => {
                        crate::report::tabular::DeflineIds::of_record(record.shown_id())
                    }
                    _ => Vec::new(),
                };
                // No id: `GetSeqIdList` throws in `SetSubjectId` (`FastaRecord::
                // subject_tabular_id`); the Seq-align id is the record's local id.
                let failed = match record {
                    Some(record) if found == Some(None) => {
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
                    id: Cow::Borrowed(id),
                    accession: Cow::Borrowed(id),
                    accession_version: Cow::Borrowed(id),
                    gi: 0,
                    all_ids,
                    titles: None,
                    length,
                    taxonomy: Cow::Owned(TabularTaxonomy::default()),
                    failed,
                    index: s_idx,
                })
            }
            TabularSubjects::Database {
                ids,
                ordinals,
                oids,
                records,
                parsed,
                sources,
                lengths,
                taxonomy,
                loads,
            } => {
                if let (Some(exception), Some(subject)) = (
                    loads.thrown(s_idx),
                    parsed.get(s_idx).and_then(Option::as_deref),
                ) {
                    return Ok(SubjectValues::thrown(subject, exception, s_idx));
                }
                let length = lengths.get(s_idx).copied().unwrap_or_default();
                let taxonomy = match taxonomy.get(s_idx) {
                    Some(taxonomy) => Cow::Borrowed(taxonomy),
                    None => Cow::Owned(TabularTaxonomy::default()),
                };
                match parsed.get(s_idx).and_then(Option::as_deref) {
                    // A sequence with parsed Seq-ids: the ids of its Bioseq and of each
                    // definition line, and the titles of its definition lines.
                    Some(subject) => {
                        let refused = |problem: &str| {
                            crate::algorithm::common::database::parsed_refusal(
                                oids.oid(ordinals.get(s_idx).copied().unwrap_or_default()),
                                &format!("has {problem}"),
                                program,
                            )
                        };
                        let ids = subject.tabular_ids(false).map_err(refused)?;
                        let all = if needs.all_deflines {
                            subject.all_tabular_ids(false).map_err(refused)?
                        } else {
                            Default::default()
                        };
                        let failed = FailedValues::of_parsed(subject, ids.is_none(), all.failed);
                        let [id, accession, accession_version] = ids.unwrap_or_default();
                        Ok(SubjectValues {
                            id: Cow::Owned(id),
                            accession: Cow::Owned(accession),
                            accession_version: Cow::Owned(accession_version),
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
                            index: s_idx,
                        })
                    }
                    // A sequence without parsed Seq-ids: the first word of its stored title (an
                    // empty one made from the name of its source; `gnl|BL_ORD_ID|n` replaced,
                    // `CShowBlastDefline::GetSeqIdList`) for every id, or where it gives none the
                    // general id (`TabularSubjectId`), and its one definition line's title.
                    None => {
                        let found = ids.get(s_idx).map(Option::as_deref);
                        let accession: &[u8] =
                            found.map_or(&b"unknown"[..], Option::unwrap_or_default);
                        let title = records
                            .get(s_idx)
                            .map_or(&[][..], |record| record.title.as_slice());
                        let oid = oids.oid(ordinals.get(s_idx).copied().unwrap_or_default());
                        // No id: `GetSeqIdList` throws for the general id
                        // (`TabularSubjectId::Unknown`), in `SetSubjectId` and `x_SetSubjectIds`.
                        let failed = (found == Some(None)).then(|| {
                            let [full, accession] =
                                crate::algorithm::common::database::failed_general_lookup(oid);
                            FailedValues {
                                align_id: [Cow::Owned(full), Cow::Owned(accession)],
                                subject_id: true,
                                all_ids: true,
                                thrown: None,
                                empty_handle: false,
                            }
                        });
                        Ok(SubjectValues {
                            id: crate::algorithm::common::database::tabular_seqid(
                                accession,
                                oid,
                                title,
                                sources
                                    .get(s_idx)
                                    .unwrap_or(&crate::report::defline::UNNAMED_SOURCE),
                                false,
                            ),
                            accession: Cow::Borrowed(accession),
                            accession_version: Cow::Borrowed(accession),
                            gi: 0,
                            all_ids: Vec::new(),
                            titles: Some(vec![title]),
                            length,
                            taxonomy,
                            failed,
                            index: s_idx,
                        })
                    }
                }
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:844-851
/// The taxonomy of the tabular fields (`staxid`, `staxids`;
/// `report::tabular::TabularTaxonomy::of_deflines`) of the reported sequences `ordinals` of `db`,
/// by rank: from the definition lines of a sequence with parsed Seq-ids (`parsed[rank]`), else
/// from its stored ones (`SeqDb::defline_set`, which keeps whether a line has a `taxid`), for
/// the ranks `ranks` (the subjects of the report's rows) that a report loads (`loads`, by rank);
/// the others are left empty: NCBI reads no definition line of a subject whose load throws (its
/// rows keep the taxonomy of the rows before) nor of one that no report reaches, and looks none
/// of their taxids up (review 1 of line DBID pfoff-2, H1 and M1). A header that does not decode
/// is refused before the report (`program` names the program).
pub(crate) fn tabular_taxonomy(
    db: &SeqDb,
    ordinals: &[u32],
    parsed: &[Option<Arc<ParsedSubject>>],
    ranks: impl IntoIterator<Item = usize>,
    loads: &crate::report::db_subject::HitLoads,
    program: &str,
) -> Result<Vec<TabularTaxonomy>> {
    let mut taxonomy = vec![TabularTaxonomy::default(); ordinals.len()];
    let mut read = vec![false; ordinals.len()];
    for rank in ranks {
        if rank >= ordinals.len() || !loads.loaded(rank) || std::mem::replace(&mut read[rank], true)
        {
            continue;
        }
        taxonomy[rank] = match parsed.get(rank).and_then(Option::as_deref) {
            Some(subject) => TabularTaxonomy::of_deflines(subject.deflines()),
            None => {
                let set = db.defline_set(ordinals[rank]).map_err(|error| {
                    crate::algorithm::common::database::database_error(error, program)
                })?;
                TabularTaxonomy::of_deflines(&set.deflines)
            }
        };
    }
    Ok(taxonomy)
}

/// Writes outfmt 20's line of the field names (`PrintProlog`, before NCBI reads the first query
/// batch), after the taxonomy database's warning of the formatter that NCBI makes for it when the
/// field list checks it and `taxdb.bti` is not found (`x_CheckTaxDB`); nothing for the other
/// kinds. `program` is the program's name in the warning (lower case), `program_version` the
/// report's (`TBLASTN 2.17.0+`). A run that finds the taxonomy database writes no warning, so
/// LOSAT's refusals (made before) come before the line.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
pub(crate) fn write_tabular_prolog(
    spec: &TabularSpec,
    taxdb_found: bool,
    program: &str,
    program_version: &str,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.kind != TabularKind::CsvWithHeader {
        return Ok(());
    }
    if spec.checks_taxdb() && !taxdb_found {
        writer.flush()?;
        diagnostics.write_all(taxdb_warning_line(program).as_bytes())?;
        diagnostics.flush()?;
    }
    TabularWriter::new(spec, program_version, None).prolog(writer)?;
    Ok(())
}

/// Writes the organism report (outfmt 18, `report::tax_report`) of the queries `queries` (title,
/// record length and HSPs in report order) of a program whose adapter has the shared subjects
/// `subjects` (`records`: a subject set's records), each query's warnings before its report. The
/// report is made before the batch's first byte (`tax_report::report_batch`), so that LOSAT's
/// refusals come before it; `names` carries the run's taxonomy database (`None`: no `taxdb.bti`
/// on the search path). `protein`: the subjects are proteins. `html`: `-html`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_tax_report(
    subjects: &TabularSubjects<'_>,
    records: &[FastaRecord],
    protein: bool,
    names: Option<&TabularNames>,
    program: &str,
    queries: &[crate::report::tax_report::TaxQuery<'_>],
    writer: &mut crate::cli::NcbiOstream<'_>,
    warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
    html: bool,
) -> Result<()> {
    use crate::report::tax_report::{report_batch, write_report};
    let report = report_batch(
        queries,
        subjects.parsed(),
        |s_idx| subjects.tax_plain(s_idx, records),
        protein,
        subjects.loads(),
        names,
        program,
        html,
    )?;
    write_report(writer, warnings, queries, &report)?;
    Ok(())
}
