//! BLASTX formatting from real final native runtime results.
use super::{
    args::ResolvedOptions, database::ReportedSubjects, input::FastaRecord, results::Hsp,
    runtime::BatchResults,
};
use crate::algorithm::tblastx::report::display_codon;
use crate::{common::GapEditOp, report::query_warnings::EMPTY_DATABASE_MESSAGE};
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
use super::empty_rows::SubjectLoc;
use crate::report::query_anchored::MaskSeek;
use crate::utils::matrix::protein_display_score;
use anyhow::{ensure, Result};
use std::{io::Write, path::Path};
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-987
#[derive(Debug)]
pub(crate) struct DisplayAlignment {
    pub qseq: String,
    pub sseq: String,
    pub gaps: usize,
    pub gap_opens: usize,
    pub qstart: i32,
    pub qend: i32,
    pub sstart: i32,
    pub send: i32,
    pub frame: i8,
    /// The columns (0-based) whose subject letter is a stored gap character (`-`): residues of
    /// the subject that the rows show as the letter of a gap (`PairwiseHit::subject_stored_gaps`).
    pub subject_stored_gaps: Vec<u32>,
}

/// The identities and positives of a display alignment's rows (`counted_alignment`): outfmt 0's hits
/// (`pairwise_hits`) count them; the tabular writer counts its rows' own
/// (`report::tabular::TabularWriter::row`), so its `alignment` makes none.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RowCounts {
    pub identity: usize,
    pub positive: usize,
}

impl RowCounts {
    // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:476-480
    /// The mismatches of `a` as NCBI prints them: fewer than none when stored gap characters of
    /// the subject face gaps of the query (equal letters, so identities, and gaps).
    pub(crate) fn mismatches(self, a: &DisplayAlignment) -> i64 {
        a.qseq.len() as i64 - self.identity as i64 - a.gaps as i64
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1348-1361
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1034-1043
pub(crate) fn query_endpoints(length: usize, frame: i8, start: i32, end: i32) -> (i32, i32) {
    let f = i32::from(frame);
    if f > 0 {
        (3 * start + f, 3 * end + f - 1)
    } else {
        (
            length as i32 - 3 * start + f + 1,
            length as i32 - 3 * end + f + 2,
        )
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:672-674
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1484-1490
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1606-1626
// This is Seq-align construction, after runtime search/pruning and e-value
pub(super) fn has_seqalign(hsp: &Hsp, gapped: bool) -> bool {
    !gapped || hsp.hsp.score != 0
}

/// Whether the hit list `list` has an HSP that the reports show (`has_seqalign`), so that
/// they read its subject's title.
pub(super) fn reported(list: &super::results::HspList, gapped: bool) -> bool {
    list.hsps.iter().any(|hsp| has_seqalign(hsp, gapped))
}

/// RPSTBLASTN: the scopes of the batches (`AnchoredBatchScopes`) where a profile of the database
/// can have the generated id (`lcl|Query_N`) of a query of the batch that reports it, which NCBI's
/// scope finds as that query (PS15a-1, PS15a-2, followup `rps_lcl_query_title.md`); `None` for
/// every other program, a subject set and queries read with `-parse_deflines`.
pub(super) fn rps_generated_scopes(
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
) -> Option<crate::blastinput::query_batch::AnchoredBatchScopes> {
    if options.rpstblastn.is_none() || !matches!(shown, ShownSubjects::Database { .. }) {
        return None;
    }
    let ids: Vec<(
        &str,
        Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
    )> = queries
        .iter()
        .map(|query| (query.internal_id.as_str(), query.believed()))
        .collect();
    let scopes = crate::blastinput::query_batch::AnchoredBatchScopes::of_ids(
        &ids,
        batches
            .iter()
            .map(|batch| batch.query_ordinal..batch.query_ordinal + batch.queries.len())
            .collect(),
    );
    scopes.generated().then_some(scopes)
}

/// The refusal of a reported list of the query `ordinal` (`lists`, those that `reported` keeps
/// and where the format differs for its profile, `differs`) whose profile has the generated id of
/// a query of its batch (`rps_generated_scopes`, `common::database::refuse_generated_query_id`),
/// with what NCBI does with it in the format asked for (`effect`).
pub(super) fn refuse_rps_generated_ids(
    scopes: &crate::blastinput::query_batch::AnchoredBatchScopes,
    shown: &ShownSubjects<'_>,
    lists: &[super::results::HspList],
    ordinal: usize,
    options: &ResolvedOptions,
    differs: impl Fn(&super::results::HspList, &crate::report::db_subject::ParsedSubject) -> bool,
    effect: &str,
) -> Result<()> {
    let ShownSubjects::Database {
        reported: names, ..
    } = shown
    else {
        return Ok(());
    };
    for list in lists.iter().filter(|list| reported(list, options.gapped)) {
        if let Some(subject) = names
            .parsed_subject(list.oid)
            .filter(|subject| differs(list, subject))
        {
            crate::algorithm::common::database::refuse_generated_query_id(
                subject,
                names.oid(list.oid),
                scopes.of_query(ordinal),
                options.program().upper,
                effect,
            )?;
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-904 (cited at
// `report::tabular::TabularSpec::subject_fetch`)
/// Whether a report whose batch's prefetch is off reads the subject of the hit list `list`: one
/// of the HSPs that it shows (`has_seqalign`) has a row that reads it (`rows`, by the engine's
/// identities, the Seq-align's `num_ident`); for `SubjectFetch::Every`, `reported`.
pub(super) fn reads_list(
    list: &super::results::HspList,
    gapped: bool,
    rows: crate::report::tabular::SubjectFetch,
) -> bool {
    list.hsps
        .iter()
        .any(|hsp| has_seqalign(hsp, gapped) && rows.reads(hsp.identity))
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:162-176
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:903-918

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:116-125
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:125-133
pub(super) fn display_query_residue(
    sequence: &[u8],
    frame: i8,
    offset: usize,
    code: &crate::utils::genetic_code::GeneticCode,
) -> u8 {
    let first = 3 * offset + frame.unsigned_abs() as usize - 1;
    let masks = std::array::from_fn(|i| {
        let position = if frame > 0 {
            first + i
        } else {
            sequence.len() - first - i - 1
        };
        let mask = b"-ACMGRSVTWYHKDBN"
            .iter()
            .position(|&c| c == sequence[position])
            .expect("validated IUPAC query residue") as u8;
        if frame > 0 {
            mask
        } else {
            ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3)
        }
    });
    display_codon(masks, code)
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:200-244
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:996-1028
pub(crate) fn alignment(
    hsp: &Hsp,
    subject: &FastaRecord,
    query_record: &FastaRecord,
    code: &crate::utils::genetic_code::GeneticCode,
) -> Result<DisplayAlignment> {
    Ok(display_alignment(hsp, subject, query_record, code, None)?.0)
}

/// `alignment` and the counts of its rows (`RowCounts`), the positives by `matrix` (the
/// formatter's matrix of the search): outfmt 0's hits (`pairwise_hits`).
pub(crate) fn counted_alignment(
    hsp: &Hsp,
    subject: &FastaRecord,
    query_record: &FastaRecord,
    code: &crate::utils::genetic_code::GeneticCode,
    matrix: crate::config::ScoringMatrix,
) -> Result<(DisplayAlignment, RowCounts)> {
    display_alignment(hsp, subject, query_record, code, Some(matrix))
}

/// The display alignment of `alignment`, with the counts of its rows where `counts` gives the
/// matrix (`counted_alignment`; none, zeros that `alignment` drops, otherwise). One function out
/// of line for both, as before the counts were split off (no NCBI counterpart: code placement
/// only).
#[inline(never)]
fn display_alignment(
    hsp: &Hsp,
    subject: &FastaRecord,
    query_record: &FastaRecord,
    code: &crate::utils::genetic_code::GeneticCode,
    counts: Option<crate::config::ScoringMatrix>,
) -> Result<(DisplayAlignment, RowCounts)> {
    let raw = &hsp.hsp;
    let (q_end, s_end) = hsp.seqalign_ends();
    let query_length = query_record.sequence.len();
    let mut qi = usize::try_from(raw.q_start)?;
    let mut si = usize::try_from(raw.s_start)?;
    let mut qseq = Vec::new();
    let mut sseq = Vec::new();
    let mut gap_opens = 0;
    // The gaps of the Dense-seg (`GetAlignLengths`): the columns of a gap op, not the `-` letters
    // of the rows (a stored gap character of the subject is one).
    let mut gaps = 0;
    let mut subject_stored_gaps = Vec::new();
    for op in display_script(hsp)? {
        let (count, consume_q, consume_s) = match op {
            GapEditOp::Sub(n) => (n, true, true),
            GapEditOp::Ins(n) => (n, true, false),
            GapEditOp::Del(n) => (n, false, true),
        };
        if !consume_q || !consume_s {
            gap_opens += 1;
            gaps += count as usize;
        }
        for _ in 0..count {
            // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:162-176
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
            let qc = if consume_q {
                display_query_residue(&query_record.sequence, raw.frame, qi, code)
            } else {
                b'-'
            };
            let sc = if consume_s {
                subject.sequence[si]
            } else {
                b'-'
            };
            if consume_s && sc == b'-' {
                subject_stored_gaps.push(u32::try_from(sseq.len())?);
            }
            qseq.push(qc);
            sseq.push(sc);
            qi += usize::from(consume_q);
            si += usize::from(consume_s);
        }
    }
    ensure!(
        qi == q_end as usize && si == s_end as usize,
        "BLASTX edit script does not reach final endpoints"
    );
    let (qstart, qend) = query_endpoints(query_length, raw.frame, raw.q_start, q_end);
    // A record cut to `-query_loc`: the reports read the frame from the record
    // (NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:411-422,
    // `seq_range::record_frame`).
    let query_place = query_record.placement();
    let frame = if query_place.offset == 0 && query_place.length == query_length {
        raw.frame
    } else {
        crate::blastinput::seq_range::record_frame(
            raw.frame > 0,
            (qstart - 1) as usize + query_place.offset,
            query_place.length,
        ) as i8
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1003-1022 (its
    // positives cited at `report::tabular::TabularWriter::new`): the tabular writer counts its
    // rows' identities and positives itself, so only outfmt 0's hits ask for them (`counts`).
    let mut row_counts = RowCounts {
        identity: 0,
        positive: 0,
    };
    if let Some(matrix) = counts {
        for (&qc, &sc) in qseq.iter().zip(&sseq) {
            if qc == sc {
                row_counts.identity += 1;
                row_counts.positive += 1;
            } else if qc != b'-' && sc != b'-' && protein_display_score(matrix, qc, sc) > 0 {
                row_counts.positive += 1;
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
    // The coordinates of records cut to `-query_loc` / `-subject_loc` move into the records by the
    // starts of their intervals (the subject's alike, `s_RemapToSubjectLoc`).
    let query_offset = query_place.offset as i32;
    let subject_offset = subject.placement().offset as i32;
    Ok((
        DisplayAlignment {
            qseq: String::from_utf8(qseq)?,
            sseq: String::from_utf8(sseq)?,
            gaps,
            gap_opens,
            qstart: qstart + query_offset,
            qend: qend + query_offset,
            sstart: raw.s_start + 1 + subject_offset,
            send: s_end + subject_offset,
            frame,
            subject_stored_gaps,
        },
        row_counts,
    ))
}

/// The operations of the display alignment of `hsp` (`alignment`, a column for each letter of an
/// operation): its edit script, or one substitution of its query span without one.
fn display_script(hsp: &Hsp) -> Result<Vec<GapEditOp>> {
    let raw = &hsp.hsp;
    Ok(if hsp.edit_script.is_empty() {
        vec![GapEditOp::Sub(u32::try_from(raw.q_end - raw.q_start)?)]
    } else {
        hsp.edit_script.clone()
    })
}

/// The values of the outfmt 0 row of `hsp`, of a query of `query_length` letters
/// (`pairwise_hits`), that the report's sorts read (`HspScores::of`), without making the row, for
/// the model of the reads (`ReportedSubjects::fetch_by_query`): the row's columns (those of
/// `display_script`, `alignment`'s rows), `alignment`'s coordinates (`query_endpoints` and the
/// subject's), the engine's identities (the row's `seqalign_num_ident`). The records' offsets of
/// a ranged query move every HSP of the query alike (a database subject has none), which no sort
/// reads, and are left out.
pub(crate) fn hsp_scores(
    hsp: &Hsp,
    query_length: usize,
) -> Result<crate::report::hit_sort::HspScores> {
    let raw = &hsp.hsp;
    let (q_end, s_end) = hsp.seqalign_ends();
    let length = display_script(hsp)?
        .iter()
        .map(|op| op.num() as usize)
        .sum();
    let (qstart, qend) = query_endpoints(query_length, raw.frame, raw.q_start, q_end);
    Ok(crate::report::hit_sort::HspScores::of_values(
        hsp.bit_score,
        hsp.evalue,
        hsp.identity,
        length,
        (usize::try_from(qstart)?, usize::try_from(qend)?),
        (usize::try_from(raw.s_start + 1)?, usize::try_from(s_end)?),
    ))
}

/// What the reports show of the subjects: `-subject`'s path (`User specified sequence
/// set`), or `-db`'s name, its `Database:` lines and the reported sequences of the query
/// batch, by rank (`database::ReportedSubjects`: ids, stored titles).
#[derive(Clone, Copy)]
pub(crate) enum ShownSubjects<'a> {
    Fasta {
        path: &'a Path,
        /// The records whose ids `-parse_deflines` read (empty without the option).
        parsed: &'a RecordsParsed,
    },
    Database {
        name: &'a str,
        pairwise: &'a crate::report::pairwise::PairwiseDatabase,
        reported: &'a ReportedSubjects,
    },
}

impl ShownSubjects<'_> {
    /// The taxonomy database's messages of the query `ordinal`'s own fetches, which follow its
    /// warnings where `PRE_FETCH_SEQS_LIMIT` turned its batch's prefetch off
    /// (`ReportedSubjects::fetch_by_query`); empty otherwise.
    pub(crate) fn query_fetches(&self, ordinal: usize) -> &[u8] {
        match self {
            ShownSubjects::Database { reported, .. } => reported.query_fetches(ordinal),
            ShownSubjects::Fasta { .. } => &[],
        }
    }

    /// The texts of the subjects that the query `ordinal`'s report reads one by one, which its
    /// writer places where it reads each (`ReportedSubjects::query_reads`); empty otherwise.
    pub(crate) fn query_reads(&self, ordinal: usize) -> &crate::report::query_warnings::QueryReads {
        static NONE: crate::report::query_warnings::QueryReads = Vec::new();
        match self {
            ShownSubjects::Database { reported, .. } => reported.query_reads(ordinal),
            ShownSubjects::Fasta { .. } => &NONE,
        }
    }
}

/// The records of a subject set whose ids `-parse_deflines` read, by subject index: as subjects
/// (`db_subject::parsed_records`), with what the tabular rows show of each (`None`: shown as a
/// record without the option).
#[derive(Default)]
pub(crate) struct RecordsParsed {
    pub(crate) subjects: Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>>,
    pub(crate) tabular: Vec<Option<RecordTabularIds>>,
}

/// `sseqid`, `sacc` and `saccver` (`None` where the lookup throws, `ParsedSubject::tabular_ids`),
/// and the `sall*` of a subject in `RecordsParsed`.
pub(crate) struct RecordTabularIds {
    /// `None` where the lookup throws (`ParsedSubject::tabular_ids`).
    pub(crate) ids: Option<[String; 3]>,
    pub(crate) all_ids: Vec<crate::report::tabular::DeflineIds>,
}

impl RecordsParsed {
    /// The subject of record `index` as outfmt 0 shows it, `None` for one shown as a record read
    /// without `-parse_deflines` (`ParsedSubject::shown_as_record`).
    pub(crate) fn shown(&self, index: usize) -> Option<&crate::report::db_subject::ParsedSubject> {
        self.subjects
            .get(index)
            .and_then(Option::as_deref)
            .filter(|subject| !subject.shown_as_record())
    }

    /// The records of `records`. `pairwise` gets the parsed subjects first, for LOSAT's refusals of
    /// them in outfmt 0 to 4 (`db_subject::check_pairwise_records`, each such format's), then come
    /// the refusals of their tabular ids.
    pub(crate) fn of(
        records: &[FastaRecord],
        show_gis: bool,
        pairwise: impl FnOnce(
            &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
        ) -> Result<()>,
    ) -> Result<Self> {
        let subjects = crate::report::db_subject::parsed_records(
            records
                .iter()
                .map(|record| (record.believed(), record.title.as_slice())),
            show_gis,
        );
        pairwise(&subjects)?;
        let refused = |subject: &crate::report::db_subject::ParsedSubject, problem: &str| {
            anyhow::anyhow!(
                "the subject record {} has {problem}, which is not supported by LOSAT's BLASTX",
                String::from_utf8_lossy(&subject.alignment_fasta())
            )
        };
        let tabular = subjects
            .iter()
            .map(|subject| {
                let Some(subject) = subject else {
                    return Ok(None);
                };
                let text = |bytes: Vec<u8>| String::from_utf8_lossy(&bytes).into_owned();
                let ids = subject
                    .tabular_ids(true)
                    .map_err(|problem| refused(subject, problem))?;
                // A record's definition line is not looked up (`AllTabularIds::failed` is false).
                Ok(Some(RecordTabularIds {
                    ids: ids.map(|ids| ids.map(text)),
                    all_ids: subject
                        .all_tabular_ids(true)
                        .map_err(|problem| refused(subject, problem))?
                        .ids,
                }))
            })
            .collect::<Result<_>>()?;
        Ok(Self { subjects, tabular })
    }
}

impl ShownSubjects<'_> {
    /// Which lookups of subject `oid`'s tabular ids fail (`report::tabular::LookupFailure`) for a
    /// report of `spec` (`None`: the default fields): a record's (with parsed ids, no
    /// `RecordTabularIds::ids`; else `subject_tabular_id`), or a database sequence's (no
    /// `ReportedSubjects::tabular_ids`; with parsed Seq-ids and a report of `sall*`,
    /// `ParsedSubject::lookup_failure`).
    pub(super) fn lookup_failure(
        &self,
        spec: Option<&crate::report::tabular::TabularSpec>,
        subject: &FastaRecord,
        oid: usize,
    ) -> Option<crate::report::tabular::LookupFailure> {
        use crate::report::tabular::LookupFailure;
        match self {
            ShownSubjects::Fasta { parsed, .. } => {
                let failed = match parsed.tabular.get(oid).and_then(Option::as_ref) {
                    Some(ids) => ids.ids.is_none(),
                    None => subject.subject_tabular_id().is_none(),
                };
                failed.then(|| LookupFailure::of_plain(false))
            }
            // A subject whose load throws where the query's batch has its prefetch off
            // (`report::db_subject::HitLoads::thrown`) fails there, before every other lookup.
            ShownSubjects::Database { reported, .. } if reported.loads().thrown(oid).is_some() => {
                Some(LookupFailure::thrown())
            }
            ShownSubjects::Database { reported, .. } => match reported.parsed_subject(oid) {
                Some(parsed)
                    if spec.is_some_and(crate::algorithm::common::tabular::reports_all_ids) =>
                {
                    parsed.lookup_failure(true).ok().flatten()
                }
                Some(_) => reported
                    .tabular_ids(oid)
                    .is_none()
                    .then_some(LookupFailure {
                        subject_id: true,
                        all_ids: false,
                        bioseq: false,
                        lineless: false,
                    }),
                None => reported
                    .tabular_ids(oid)
                    .is_none()
                    .then(|| LookupFailure::of_plain(true)),
            },
        }
    }

    /// The `-db` database's `Database:` lines.
    fn pairwise_database(&self) -> Option<&crate::report::pairwise::PairwiseDatabase> {
        match self {
            ShownSubjects::Fasta { .. } => None,
            ShownSubjects::Database { pairwise, .. } => Some(pairwise),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:289-290
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
/// The reports of one query batch, the `queries_before` queries of the batches before it written
/// already. `stops`: the run stops after this batch (a later batch fails before its search), so
/// a document that holds each report back never writes the last query's
/// (`blast_xml::written_queries`). `xml_matrix`: where the run's XML formats fail for a matrix
/// name NCBI's XML does not know (`blast_xml::xml_matrix_failure`, decided once for the run, at the
/// run's query of `QueryScopes::xml_failing_query`).
/// `pool`: the run's pool where no search runs on it while the batch is written
/// (`native::BatchRender`), whose threads format the pairwise reports (`x_pairwise_par::write`).
pub(crate) fn render_queries(
    writer: &mut crate::cli::NcbiOstream<'_>,
    stderr: &mut dyn Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &mut [BatchResults],
    queries_before: usize,
    stops: bool,
    xml2: Option<&mut super::blast_xml::RunXml2>,
    xml_matrix: Option<crate::blastinput::query_batch::XmlFailingQuery>,
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    sort_lists(batches, SubjectLoc::of(shown));
    render_queries_with(
        writer,
        stderr,
        options,
        queries,
        subjects,
        shown,
        batches,
        queries_before,
        stops,
        xml2,
        xml_matrix,
        pool,
        None,
    )
}

/// What LOSAT Web's `run_local` (`local.rs`) adds to the reports of one batch in one format.
pub(super) struct LocalReport<'a, 'p> {
    /// The format's observer (`api::local_blast::FormatProbe`).
    pub probe: Option<&'a mut crate::api::local_blast::FormatProbe<'p>>,
    /// The index of the batch's first reported HSP in the run's hit list.
    pub first_hsp: crate::api::local_blast::HspIndex,
    /// The batch's hits (`pairwise_hits`), when they are made already.
    pub hits: Option<&'a [crate::report::pairwise::PairwiseHit]>,
    /// What a refusal of the format's checks does (`FormatRefusal::check`); without a local
    /// report it fails the call.
    pub refusal: &'a mut crate::api::local_blast::FormatRefusal,
    /// The format's place among the search's formats (`QueryWarnings::for_report`); 0 without a
    /// local report (the command line's one format).
    pub report: usize,
    /// Where the format's `-html` posts go apart from its standard error
    /// (`QueryWarnings::with_html_posts`: the diagnostics of a later format, the call's receiver
    /// of a later search's, `local::write_run`).
    pub html_posts: Option<&'a mut dyn Write>,
}

/// Puts each batch of `batches`, whose subjects are `loc`'s, in the order of every report
/// (`BatchResults::put_in_report_order`: sorted by e-value, `-subject`'s HSPs mapped, once).
/// Each batch's report calls it before any format reads the lists (`render_queries`, the
/// archive's `native::ArchiveBatch::write`, `run_local`'s batch).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
pub(super) fn sort_lists(batches: &mut [BatchResults], loc: SubjectLoc) {
    for batch in batches.iter_mut() {
        batch.put_in_report_order(loc);
    }
}

/// The number of HSPs that the reports of `batches` show (`has_seqalign`).
pub(super) fn reported_hsps(batches: &[BatchResults], gapped: bool) -> usize {
    batches
        .iter()
        .flat_map(|batch| batch.queries.iter().flatten())
        .map(|list| list.hsps.iter().filter(|h| has_seqalign(h, gapped)).count())
        .sum()
}

/// LOSAT's refusals of the tabular report (`options.tabular`, `None`: the default fields) of the
/// batches' rows whose subject lookup fails where NCBI has no result
/// (`report::tabular::check_failed_lookups`; each query of a batch is a result set), with the
/// query (its index in `queries`, the call's: the batch's `query_ordinal` and the set) whose row
/// NCBI stops at.
fn check_failed_lookups(
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
) -> std::result::Result<(), (usize, anyhow::Error)> {
    if !crate::report::tabular::may_stop_at_failed_lookup(options.tabular) {
        return Ok(());
    }
    for batch in batches {
        crate::report::tabular::check_failed_lookups(
            options.tabular,
            (0..batch.queries.len()).map(|qidx| {
                batch
                    .shown_lists(qidx, options)
                    .iter()
                    .filter(|list| {
                        list.hsps
                            .iter()
                            .any(|hsp| has_seqalign(hsp, options.gapped))
                    })
                    .map(|list| list.oid)
            }),
            |oid| shown.lookup_failure(options.tabular, &subjects[oid], oid),
            |oid, failure| match shown {
                ShownSubjects::Database { reported, .. } => {
                    crate::report::tabular::failed_database_subject(
                        reported.oid(oid),
                        reported.loads().thrown(oid),
                        failure,
                    )
                }
                ShownSubjects::Fasta { .. } => {
                    crate::report::tabular::failed_record_title(oid, subjects[oid].shown_id())
                }
            },
            options.program().upper,
            options.taxdb_found,
        )
        .map_err(|error| {
            let set = error
                .downcast_ref::<crate::report::tabular::FailedLookupRefusal>()
                .map_or(0, crate::report::tabular::FailedLookupRefusal::set);
            (batch.query_ordinal + set, error)
        })?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
/// What NCBI writes to `stderr` before it stops at the tabular row of the query `stop` (its index
/// in `queries`) that LOSAT refuses (`error`): the batch's warning of `-subject`'s prefetch, then
/// for each query up to `stop` what the tabular report writes to the standard error of its report
/// (`tabular::write_blastx_tabular`: its warnings and fetch texts, its formatter's warning, its
/// subject reads; `FailedLookupRefusal::posts_before_stop`, whose reads name a subject by its
/// `oid`). The batch's reader messages and `-db` prefetch texts came before the report.
#[allow(clippy::too_many_arguments)]
fn posts_before_lookup_stop(
    stderr: &mut dyn Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
    stop: usize,
    error: &anyhow::Error,
) -> Result<()> {
    let Some(refusal) = error.downcast_ref::<crate::report::tabular::FailedLookupRefusal>() else {
        return Ok(());
    };
    if !options.gapped
        && SubjectLoc::of(shown) == SubjectLoc::Interval
        && super::native::batch_prefetch(batches, options, SubjectLoc::Interval).fails
    {
        stderr.write_all(&crate::algorithm::common::database::prefetch_warning(
            options.program().lower,
        ))?;
    }
    let (mut posts, mut reads) = (Vec::new(), Vec::new());
    for batch in batches {
        for qidx in 0..batch.queries.len() {
            let ordinal = batch.query_ordinal + qidx;
            if ordinal > stop {
                break;
            }
            let mut text = query_warning(options, &queries[ordinal], batch.search_skipped);
            text.extend_from_slice(shown.query_fetches(ordinal));
            posts.push(crate::report::diag_log::Posts::of_lines(text));
            reads.push(shown.query_reads(ordinal).clone());
        }
    }
    refusal
        .posts_before_stop(&posts, &reads, |oid| oid as u32)
        .write(stderr, None)?;
    Ok(())
}

/// `render_queries`, with what `run_local` adds (`local`), of `batches` whose lists are in the
/// order of the reports (`sort_lists`).
#[allow(clippy::too_many_arguments)]
pub(super) fn render_queries_with(
    writer: &mut crate::cli::NcbiOstream<'_>,
    stderr: &mut dyn Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &mut [BatchResults],
    queries_before: usize,
    stops: bool,
    xml2: Option<&mut super::blast_xml::RunXml2>,
    xml_matrix: Option<crate::blastinput::query_batch::XmlFailingQuery>,
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
    local: Option<LocalReport<'_, '_>>,
) -> Result<()> {
    let (probe, first_hsp, made_hits, mut refusal, html_posts) = match local {
        Some(local) => (
            local.probe,
            local.first_hsp,
            local.hits,
            Some(local.refusal),
            (local.report, local.html_posts),
        ),
        None => (None, 0, None, None, (0, None)),
    };
    // A check of the format's report before its rows: `false` where the format fails alone
    // (`FormatRefusal::check`), and writes nothing of the batch.
    let mut writes = |checked: Result<()>| match refusal.as_deref_mut() {
        Some(refusal) => refusal.check(checked),
        None => checked.map(|()| true),
    };
    if !options.gapped && options.outfmt != 18 && SubjectLoc::of(shown) == SubjectLoc::Interval {
        // The tabular reports (`ResolvedOptions::tabular_report`) write a subject's only HSP of
        // length 0 as NCBI's formatter reads its mapped alignment.
        let tabular = options
            .tabular_report()
            .is_some()
            .then_some(super::empty_rows::MappedTables { queries, subjects });
        if !writes(super::empty_rows::refuse_mapped_without_scores(
            batches,
            options.gapped,
            tabular,
        ))? {
            return Ok(());
        }
    }
    // The rows of subjects whose lookup fails where NCBI has no result
    // (`report::tabular::check_failed_lookups`), before the batch's tabular report and before
    // the format writes a warning of the batch (the prefetch's below, the queries'), so that a
    // format that fails alone leaves them to the next (`local::write_run`); a refusal that fails
    // the call comes after those warnings up to the row's query (`posts_before_lookup_stop`).
    if options
        .tabular_report()
        .is_some_and(crate::report::tabular::TabularSpec::writes_rows)
    {
        let (stop, checked) = match check_failed_lookups(options, subjects, shown, batches) {
            Ok(()) => (None, Ok(())),
            Err((stop, error)) => (Some(stop), Err(error)),
        };
        match writes(checked) {
            Ok(true) => {}
            Ok(false) => return Ok(()),
            Err(error) => {
                if let Some(stop) = stop {
                    posts_before_lookup_stop(
                        stderr, options, queries, shown, batches, stop, &error,
                    )?;
                }
                return Err(error);
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:286-287 (cited at
    // `native::reported_batch`)
    // `-subject`'s prefetch before the batch's reports, as `-db`'s (`native::batch_prefetch`): a
    // subject's Seq-align whose empty rows or remaining leading subject gap the mapper leaves
    // makes it write its warning.
    if SubjectLoc::of(shown) == SubjectLoc::Interval
        && super::native::batch_prefetch(batches, options, SubjectLoc::Interval).fails
    {
        stderr.write_all(&crate::algorithm::common::database::prefetch_warning(
            options.program().lower,
        ))?;
    }
    if options.pairwise_view().is_some() {
        return write_pairwise(
            writer, stderr, options, queries, subjects, shown, batches, made_hits, probe,
            first_hsp, pool, html_posts,
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1432-1444
    // BLAST XML (outfmt 5): the Iterations of the batch's queries; the header comes with the
    // first query of the run and the end of the document with the epilog. The XML2 and JSON
    // documents (outfmt 15 and 16) are the reports of the batch's queries through the run's
    // writer, whose start comes with the prolog and whose end with the epilog.
    if options.outfmt == 5 {
        return super::blast_xml::write_blastx_xml(
            writer,
            stderr,
            &mut crate::report::blast_xml::BlastXmlWriter::new(queries_before),
            options,
            queries,
            subjects,
            shown,
            batches,
            queries_before,
            stops,
            xml_matrix,
        );
    }
    if let Some(xml2) = xml2 {
        return super::blast_xml::write_blastx_xml(
            writer,
            stderr,
            xml2,
            options,
            queries,
            subjects,
            shown,
            batches,
            queries_before,
            stops,
            xml_matrix,
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
    // The Seq-annot of each query with alignments (outfmt 8, 9 and 12), of the lists that
    // `PruneSeqalign` keeps (`BatchResults::shown_lists`); the structured report returns before
    // the query's warnings (blast_format.cpp:1431-1452, `ReportFormat::writes_query_warnings`).
    if let Some(serial) =
        crate::blastinput::app::ReportFormat::from_number(i32::from(options.outfmt))
            .and_then(crate::report::serial::SerialFormat::of)
    {
        let names = match shown {
            ShownSubjects::Fasta { .. } => super::seqalign::SeqalignSubjects::Records,
            ShownSubjects::Database {
                pairwise, reported, ..
            } => super::seqalign::SeqalignSubjects::Database {
                reported,
                title: &pairwise.definition,
            },
        };
        return super::seqalign::write_blastx_seqalign(
            serial, writer, stderr, options, queries, subjects, &names, batches,
        );
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:78-81 (cited at
    // `report::tabular::TabularSpec::std`)
    // The tabular reports: the default outfmt 6 and 7 are the `std` fields and a tab of the
    // formatter (`CBlastTabularInfo`) that writes any field list, so one writer writes them all
    // (`tabular::write_blastx_tabular`, the shared rows of `report::tabular`) with the run's report
    // (`ResolvedOptions::tabular_report`); every other format returned above.
    let spec = options.tabular_report().ok_or_else(|| {
        anyhow::anyhow!(
            "internal error: outfmt {} has no tabular report",
            options.outfmt
        )
    })?;
    // RPSTBLASTN's organism report (outfmt 18) and comma-separated values with a header line
    // (outfmt 20, item PS15a-2): a profile whose Seq-align id is the generated id of a query of
    // its batch is that query for NCBI's formatter, whose rows then differ or stop with an
    // exception (`TabularSpec::refuses_generated_id`). Refused before any byte, as the documents
    // above; outfmt 6, 7 and 10 of the same rows are line RG item X7's.
    if spec.refuses_generated_id(crate::blastinput::app::Program::Blastx, 1, true) {
        if let Some(scopes) = rps_generated_scopes(options, queries, shown, batches) {
            for batch in batches.iter() {
                for qidx in 0..batch.queries.len() {
                    refuse_rps_generated_ids(
                        &scopes,
                        shown,
                        batch.shown_lists(qidx, options),
                        batch.query_ordinal + qidx,
                        options,
                        |list, subject| {
                            let taxid = subject.first_source_taxid().is_some();
                            list.hsps.iter().any(|hsp| {
                                has_seqalign(hsp, options.gapped)
                                    && spec.refuses_generated_id(
                                        crate::blastinput::app::Program::Blastx,
                                        hsp.identity,
                                        taxid,
                                    )
                            })
                        },
                        crate::algorithm::common::database::TABLE_EFFECT,
                    )?;
                }
            }
        }
    }
    super::tabular::write_blastx_tabular(
        writer, stderr, options, spec, queries, subjects, shown, batches, probe, first_hsp,
    )
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1420
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1578
#[allow(clippy::too_many_arguments)]
fn write_pairwise(
    writer: &mut impl crate::report::ostream::ReportWrite,
    stderr: &mut dyn Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
    made_hits: Option<&[crate::report::pairwise::PairwiseHit]>,
    probe: Option<&mut crate::api::local_blast::FormatProbe<'_>>,
    first_hsp: crate::api::local_blast::HspIndex,
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
    (report_index, html_posts): (usize, Option<&mut dyn Write>),
) -> Result<()> {
    use crate::report::pairwise::{
        write_blastx_pairwise_report, BlastpPairwiseQuery, BlastxPairwiseOptions, PairwiseConfig,
    };
    let mut qinfo = Vec::new();
    // The posts before each query's report (`QueryWarnings`: after the report written so far).
    let mut warnings = Vec::new();
    // The texts of the subjects that each query's report reads one by one, which the writer
    // places where it reads each (`ShownSubjects::query_reads`).
    let mut reads = Vec::new();
    for batch in batches {
        for qidx in 0..batch.queries.len() {
            let absolute = batch.query_ordinal + qidx;
            let query = &queries[absolute];
            let first = (qidx * 6..qidx * 6 + 6).find(|&c| batch.parameters[c].valid);
            let stat = first.map(|c| &batch.parameters[c]);
            qinfo.push(BlastpPairwiseQuery {
                // The footer's Karlin blocks (`write_translated_query_footer`).
                valid: first.is_some(),
                batch_skipped: batch.search_skipped,
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:630-637
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-743
                query_name: query.query_line().into_owned(),
                query_length: query.placement().length,
                ungapped_karlin: stat.map_or(batch.parameters[qidx * 6].ungapped, |p| p.ungapped),
                effective_search_space: stat.map_or(0, |p| p.search_space),
            });
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
            let mut query_warnings = query_warning(options, query, batch.search_skipped);
            query_warnings.extend_from_slice(shown.query_fetches(batch.query_ordinal + qidx));
            warnings.push(crate::report::diag_log::Posts::of_lines(query_warnings));
            reads.push(shown.query_reads(batch.query_ordinal + qidx).clone());
        }
    }
    // Under `-html` the rows' masks of each hit (`HtmlRowMasks`, the spans of outfmt 0), made with
    // the rows: the hits made for the caller (`made_hits`, the same rows) are made again with them.
    let mut html_masks = options.html.is_some().then(Vec::new);
    let made;
    let hits = match made_hits {
        Some(hits) if html_masks.is_none() => hits,
        _ => {
            made = pairwise_hits(
                options,
                queries,
                subjects,
                shown,
                batches,
                html_masks.as_mut(),
            )?;
            &made
        }
    };
    let report = pairwise_metadata(options, subjects, shown);
    // The description line of each subject with hits: the title that `GenerateDefline` makes
    // (`Unknown` where it throws), or a sequence's with parsed Seq-ids
    // (`database_description_label`), as the other programs' (D2 of line RG T7: no longer the
    // record's first word and rest); the others are not shown.
    let ids: Vec<crate::report::pairwise::DescriptionLabel> = match shown {
        ShownSubjects::Fasta { parsed, .. } => {
            let mut reported = vec![false; subjects.len()];
            for hit in hits {
                if let Some(seen) = reported.get_mut(hit.hit.s_idx as usize) {
                    *seen = true;
                }
            }
            subjects
                .iter()
                .enumerate()
                .map(|(index, s)| {
                    if !reported[index] {
                        return Default::default();
                    }
                    crate::report::pairwise::database_description_label(
                        &s.title,
                        &crate::report::defline::UNNAMED_SOURCE,
                        parsed.shown(index),
                        true,
                        false,
                    )
                })
                .collect()
        }
        ShownSubjects::Database { reported, .. } => reported.descriptions(),
    };
    // The titles of a subject set's records, whose alignment headings `GenerateDefline` makes.
    let record_titles: Vec<&[u8]> = match shown {
        ShownSubjects::Fasta { .. } => subjects.iter().map(|s| s.title.as_slice()).collect(),
        ShownSubjects::Database { .. } => Vec::new(),
    };
    let database = match shown {
        ShownSubjects::Fasta { .. } => None,
        ShownSubjects::Database {
            pairwise, reported, ..
        } => Some((
            *pairwise,
            reported.stored_titles(),
            reported.parsed(),
            reported.sources(),
        )),
    };
    let records_parsed = match shown {
        ShownSubjects::Fasta { parsed, .. } => parsed.subjects.as_slice(),
        ShownSubjects::Database { .. } => &[],
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:247-249
    let view = options
        .pairwise_view()
        .expect("a pairwise report is outfmt 0 to 4");
    // The query-anchored formats (outfmt 1 to 4, line RG T3c), before the call's first byte: a
    // reported database sequence that NCBI's scope finds as a query of its batch, by its Seq-align
    // id (a generated `lcl|Query_N`) or, in a description line, by its display id
    // (`common::database::refuse_query_ids_of_anchored_reports`, as BLASTN's and TBLASTN's): NCBI
    // shows the query in its place.
    if let (true, ShownSubjects::Database { reported, .. }) = (
        view != crate::report::pairwise::PairwiseView::Pairwise,
        shown,
    ) {
        let ids: Vec<(
            &str,
            Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
        )> = queries
            .iter()
            .map(|query| (query.internal_id.as_str(), query.believed()))
            .collect();
        let scopes = crate::blastinput::query_batch::AnchoredBatchScopes::of_ids(
            &ids,
            batches
                .iter()
                .map(|batch| batch.query_ordinal..batch.query_ordinal + batch.queries.len())
                .collect(),
        );
        let mut reported_ranks: Vec<Vec<u32>> = vec![Vec::new(); queries.len()];
        for hit in hits {
            let ranks = &mut reported_ranks[hit.hit.q_idx as usize];
            if !ranks.contains(&hit.hit.s_idx) {
                ranks.push(hit.hit.s_idx);
            }
        }
        let parsed = reported.parsed();
        crate::algorithm::common::database::refuse_query_ids_of_anchored_reports(
            queries.len(),
            |q_idx| std::mem::take(&mut reported_ranks[q_idx]),
            |q_idx| {
                Ok(crate::report::pairwise::blastp_described_subjects(
                    view.report_sort(options.report_sort()),
                    None,
                    options.num_descriptions as usize,
                    q_idx,
                    hits,
                    |previous, s_idx| {
                        crate::report::db_subject::rank_joins_as_copy(
                            parsed,
                            previous as usize,
                            s_idx as usize,
                        )
                    },
                ))
            },
            parsed,
            |rank| reported.oid(rank as usize),
            |q_idx| scopes.of_query(q_idx),
            scopes.generated(),
            options.program().upper,
        )?;
    }
    let config = PairwiseConfig {
        program: "blastx".into(),
        protein_matrix: options.scoring,
        line_length: options.formatting.line_length(),
        // `-sorthits` for every view, `-sorthsps` for outfmt 0 only (NCBI's warning for the others,
        // `PairwiseView::report_sort`).
        sort: view.report_sort(options.report_sort()),
        view,
        html: options.html.as_ref().map(|links| {
            let oids: &[crate::seqdb::Oid] = match shown {
                ShownSubjects::Database { reported, .. } => reported.oids(),
                ShownSubjects::Fasta { .. } => &[],
            };
            crate::report::html_ids::HtmlReport::new(links, oids.iter().copied())
        }),
        ..PairwiseConfig::default()
    };
    // The sequences of the query-anchored sections of outfmt 1 to 4 (line RG T3c), made once.
    let anchored = (view != crate::report::pairwise::PairwiseView::Pairwise).then(|| {
        super::anchored::AnchoredSequences {
            queries,
            subjects,
            database: match shown {
                ShownSubjects::Fasta { .. } => None,
                ShownSubjects::Database { reported, .. } => {
                    Some((reported.parsed(), reported.oids()))
                }
            },
            show_gis: options.show_gis,
            query_masks: super::anchored::AnchoredSequences::query_masks(queries, batches),
            subject_masks: match db_subject_masks(shown, hits) {
                Some(masks) => masks,
                None if subjects.iter().any(FastaRecord::has_lowercase_runs) => {
                    super::x_soft_subject::x_subject_masks(hits, subjects)
                }
                None => Default::default(),
            },
            code: crate::algorithm::tblastx::report::DisplayCode::new(
                &crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode),
            ),
        }
    });
    // The queries' gapped K where a batch gives its own (`BatchResults::gapped_k`), in the
    // reports' order of the queries; the table's K elsewhere.
    let query_gapped_k: Option<Vec<f64>> = batches
        .iter()
        .any(|batch| batch.gapped_k.is_some())
        .then(|| {
            batches
                .iter()
                .flat_map(|batch| match &batch.gapped_k {
                    Some(gapped_k) => gapped_k.clone(),
                    None => vec![report.gapped_karlin.k; batch.queries.len()],
                })
                .collect()
        });
    // The warnings of each query and the texts of its subjects' reads, at the shared writer's places
    // (`QueryWarnings::before_query`, `subject_read`).
    let mut warnings =
        crate::report::query_warnings::QueryWarnings::of_posts(&warnings, None, stderr)
            .with_subject_reads(&reads)
            .for_report(report_index)
            .with_html_posts(html_posts.map(|posts| posts as &mut dyn Write));
    // A subject set's subjects are all loaded.
    let every_loaded = crate::report::db_subject::HitLoads::default();
    write_blastx_pairwise_report(
        hits,
        writer,
        Some(&mut warnings),
        &config,
        &qinfo,
        &ids,
        &report,
        &BlastxPairwiseOptions {
            word_threshold: options.threshold,
            gapped: options.gapped,
            sum_stats: options.sum_stats,
            num_descriptions: options.num_descriptions as usize,
            num_alignments: options.num_alignments as usize,
        },
        database,
        records_parsed,
        &record_titles,
        query_gapped_k.as_deref(),
        probe,
        first_hsp,
        anchored
            .as_ref()
            .map(|sequences| sequences as &dyn crate::report::query_anchored::ReportSequences),
        pool.and_then(crate::report::x_pairwise_par::pool),
        html_masks.as_deref(),
        match shown {
            ShownSubjects::Fasta { .. } => &every_loaded,
            ShownSubjects::Database { reported, .. } => reported.loads(),
        },
    )
    .map_err(crate::report::pairwise::empty_rows_stop)?;
    Ok(())
}

/// The hits of the HSPs that the reports of `batches` show (`has_seqalign`), in the order of the
/// reports (the lists sorted by e-value, `sort_lists`); `q_idx` is the query's index in `queries`.
/// Outfmt 0 writes them (`write_pairwise`), and LOSAT Web's `run_local` hands them to its caller.
/// Under `-html` (`html_masks`) each hit's rows' masks are added to `html_masks` in the same order
/// (`report::pairwise::HtmlRowMasks`).
pub(super) fn pairwise_hits(
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
    mut html_masks: Option<&mut Vec<crate::report::query_anchored::HspMasks>>,
) -> Result<Vec<crate::report::pairwise::PairwiseHit>> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    use crate::common::Hit;
    use crate::report::pairwise::PairwiseHit;
    let mut hits = Vec::new();
    // The rows' masks, kept from row to row (`query_row_masks`).
    let mut row_masks = crate::report::query_anchored::HspMasks::default();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let absolute = batch.query_ordinal + qidx;
            let query = &queries[absolute];
            // The query's masks (as the query-anchored reports read them), by frame.
            let query_masks = if lists.is_empty() {
                Vec::new()
            } else {
                masks_by_frame(&super::anchored::query_record_masks(batch, qidx, query))
            };
            for list in lists {
                let subject = &subjects[list.oid];
                for hsp in &list.hsps {
                    if !has_seqalign(hsp, options.gapped) {
                        continue;
                    }
                    let (a, counts) =
                        counted_alignment(hsp, subject, query, &code, options.scoring)?;
                    let raw = &hsp.hsp;
                    // The report's hit is the HSP's Seq-align, so its offsets end where the
                    // Seq-align does (`Hsp::seqalign_ends`); their one reader here is
                    // `PairwiseHit::rows_empty`, NCBI's `CSeq_align::GetSeqRange` on these rows.
                    let (q_end, s_end) = hsp.seqalign_ends();
                    let len = a.qseq.len();
                    let hit = Hit {
                        identity: 100.0 * counts.identity as f64 / len as f64,
                        length: len,
                        mismatch: counts.mismatches(&a),
                        gapopen: a.gap_opens,
                        q_start: usize::try_from(a.qstart)?,
                        q_end: usize::try_from(a.qend)?,
                        s_start: usize::try_from(a.sstart)?,
                        s_end: usize::try_from(a.send)?,
                        e_value: hsp.evalue,
                        bit_score: hsp.bit_score,
                        num_ident: counts.identity,
                        query_frame: i32::from(a.frame),
                        query_length: query.placement().length,
                        q_idx: u32::try_from(absolute)?,
                        s_idx: u32::try_from(list.oid)?,
                        raw_score: raw.score,
                        sort_query_offset: usize::try_from(raw.q_start)?,
                        sort_query_end: usize::try_from(q_end)?,
                        sort_subject_offset: usize::try_from(raw.s_start)?,
                        sort_subject_end: usize::try_from(s_end)?,
                        has_sort_offsets: true,
                        gap_info: Some(hsp.edit_script.clone()),
                        num_positives: counts.positive,
                    };
                    // The row holds its masks in lower case; the lists are left to the next
                    // row (`query_anchored::output_seq`), or copied under `-html`.
                    if html_masks.is_some() {
                        row_masks.clear();
                    }
                    let query_seq = query_row_masks(
                        a.qseq,
                        (a.qstart, a.qend),
                        a.frame,
                        &query_masks,
                        &mut row_masks,
                    );
                    if let Some(table) = html_masks.as_deref_mut() {
                        table.push(row_masks.copy_rows());
                    }
                    hits.push(PairwiseHit {
                        hit,
                        query_seq: Some(query_seq),
                        subject_seq: Some(a.sseq),
                        query_frame: Some(a.frame),
                        subject_frame: None,
                        positives: Some(counts.positive),
                        gaps: Some(a.gaps),
                        subject_length: Some(subject.placement().length),
                        // A database sequence's description line is made from its stored
                        // title alone (`ReportedSubjects::descriptions`).
                        subject_title: match shown {
                            ShownSubjects::Fasta { parsed, .. }
                                if parsed.shown(list.oid).is_some() =>
                            {
                                None
                            }
                            // The library's `Hit::subject_title` (no report reads it),
                            // `protein product` of an empty title.
                            ShownSubjects::Fasta { .. } if subject.title.is_empty() => {
                                Some("protein product".to_owned())
                            }
                            ShownSubjects::Fasta { .. } => {
                                crate::blastinput::fasta_reader::hit_subject_title(&subject.title)
                            }
                            ShownSubjects::Database { .. } => None,
                        },
                        comp_adjust_method: Some(u8::try_from(hsp.composition_method)?),
                        sum_n: Some(hsp.num),
                        subject_stored_gaps: a.subject_stored_gaps,
                        query_stored_gaps: Vec::new(),
                        query_case: Default::default(),
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:530-531
                        // The Seq-align's `num_ident` (`seqalign::query_records`): the engine's
                        // identities on the unmasked query frame (`statistics::identity_and_positive`),
                        // not those of the rows (`a.identity`), which the display's translation
                        // and the subject's stored gap characters make.
                        seqalign_num_ident: hsp.identity,
                    });
                    // The model of the reads sorts by `hsp_scores` of the batch's query length
                    // without the rows (`ReportedSubjects::fetch_by_query`): the row's values
                    // where the records' offsets move nothing.
                    debug_assert!(
                        query.placement().offset != 0
                            || subject.placement().offset != 0
                            || hsp_scores(hsp, batch.prepared.original_lengths[qidx]).ok()
                                == hits.last().map(crate::report::hit_sort::HspScores::of),
                        "BLASTX hsp_scores differs from its row"
                    );
                }
            }
        }
    }
    // `-lcase_masking`: the subject rows show their records' lower case too (x_soft_subject.rs), in
    // outfmt 0 and in the hits that LOSAT Web's `run_local` hands out, as the query rows do.
    // A `-subject_loc` subject can show runs of its record outside the interval
    // (`FastaRecord::archive_lowercase_runs`). `-db_soft_mask` and `-db_hard_mask`: the stored
    // masks by NCBI's Seq-id (`ReportedSubjects::shown_db_masks`).
    if let Some(masks) = db_subject_masks(shown, &hits) {
        super::x_soft_subject::x_set_subject_row_masks(&mut hits, &masks, html_masks);
    } else if subjects.iter().any(FastaRecord::has_lowercase_runs) {
        let kept = super::x_soft_subject::x_subject_masks(&hits, subjects);
        super::x_soft_subject::x_set_subject_row_masks(&mut hits, &kept, html_masks);
    }
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:936-951
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2500-2503
// The masks a row shows are those labelled with the row's frame in the record; a mask keeps
// the frame of the context it was found in, which counts from the ends of the searched letters.
// With `-query_loc` the two can differ, and NCBI then shows the masks of another frame's context
// on the row (the TBLASTX display, `tblastx::report::shown_query_masks`).
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:593-595
// A position's column in its segment is its distance from the segment's first position divided
// by the row's width (three nucleotides a residue here), counted from the segment's end on the
// minus strand.
/// The outfmt 0 query row `row` of a display alignment (`DisplayAlignment::qseq`, of the 1-based
/// record positions `ends`, `qstart` and `qend`, on the frame `frame`) with its masks in lower
/// case, which it sets as row 0 of `masks` (`HspMasks`, the rows' builder's, which keeps its
/// buffers from row to row): the query's masks of the row's frame (`query_masks`,
/// `masks_by_frame`; record positions) at the row's columns of `x_FillLocList`
/// (`query_anchored::fill_row_loc_list`: the row's nucleotides from its first to its last, three
/// a residue, on the frame's strand), in the masks' order, which `x_OutputSeq` puts in lower
/// case. A mask of another frame never shows on the row (`x_OutputSeq` passes over it), so the
/// row does not hold it. The row comes back so, as the hit records that the library API and
/// LOSAT Web hand out show it. An HSP without query letters (an ungapped HSP of length 0,
/// `-window_size 0`) has no mask.
fn query_row_masks(
    row: String,
    ends: (i32, i32),
    frame: i8,
    query_masks: &[(i8, MaskSeek<u32>)],
    masks: &mut crate::report::query_anchored::HspMasks,
) -> String {
    let Some((_, seek)) = query_masks
        .iter()
        .find(|(masks_frame, _)| *masks_frame == frame)
    else {
        return row;
    };
    if row.is_empty() {
        return row;
    }
    let (first, last) = (ends.0.min(ends.1), ends.0.max(ends.1));
    let range = ((first - 1) as u32, (last - 1) as u32);
    masks.add_meeting(seek, range.0, range.1, |from, to| (from, to, frame));
    let letters = row.as_bytes();
    masks.set_row_met(
        0,
        frame,
        letters,
        |column| letters[column] != b'-',
        range,
        frame > 0,
        3,
    );
    if masks.row(0).is_none() {
        return row;
    }
    let mut row = row.into_bytes();
    masks.lower_row(0, &mut row);
    String::from_utf8(row).expect("ASCII alignment row")
}

/// The masks `masks` of a query (every frame's, `anchored::query_record_masks`) by frame, each
/// frame's in their order, arranged for the scans of the rows (`MaskSeek`): a row reads only its
/// frame's masks (`query_row_masks`).
fn masks_by_frame(masks: &[crate::report::query_anchored::Mask]) -> Vec<(i8, MaskSeek<u32>)> {
    let mut frames: Vec<i8> = Vec::new();
    for &(_, _, frame) in masks {
        if !frames.contains(&frame) {
            frames.push(frame);
        }
    }
    frames
        .into_iter()
        .map(|frame| {
            let of_frame = masks
                .iter()
                .filter(|mask| mask.2 == frame)
                .map(|&(from, to, _)| (from, to));
            (frame, MaskSeek::new(of_frame))
        })
        .collect()
}

/// The subject masks of a search with `-db_soft_mask` or `-db_hard_mask` that the rows `hits`
/// show (`ReportedSubjects::shown_db_masks`); `None` for any other search.
fn db_subject_masks(
    shown: &ShownSubjects<'_>,
    hits: &[crate::report::pairwise::PairwiseHit],
) -> Option<std::collections::HashMap<(u32, u32), Vec<(i32, i32)>>> {
    match shown {
        ShownSubjects::Database { reported, .. } => reported.shown_db_masks(hits),
        ShownSubjects::Fasta { .. } => None,
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:226-244
fn pairwise_metadata(
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
) -> crate::report::pairwise::BlastpPairwiseReport {
    use crate::{
        config::ProteinScoringSpec,
        report::pairwise::BlastpPairwiseReport,
        stats::{spouge::lookup_protein_gumbel_params, tables::lookup_protein_params},
    };
    let scoring = ProteinScoringSpec {
        matrix: options.scoring,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    };
    // A database's lines come from `ShownSubjects::pairwise_database`; these fields serve
    // `-subject` (the Gumbel block's printed numbers do not depend on the length).
    let (database_name, num_sequences, total) = match shown {
        // NCBI's `Int8` counts (`PairwiseDatabase`); a FASTA set's records are held in memory.
        ShownSubjects::Fasta { path, .. } => (
            crate::report::pairwise::subject_set_name(path),
            subjects.len() as i64,
            subjects.iter().map(|s| s.sequence.len() as i64).sum(),
        ),
        ShownSubjects::Database { name, pairwise, .. } => (
            name.as_bytes().to_vec(),
            pairwise.num_sequences,
            pairwise.total_letters,
        ),
    };
    BlastpPairwiseReport {
        version: crate::wording::REPORT_VERSION.into(),
        database_name,
        database_num_sequences: num_sequences,
        database_total_letters: total,
        matrix_name: Some(options.matrix.clone()),
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
        word_threshold: 0.0,
        window_size: options.window_size,
        // An ungapped search's gap costs need not be in the matrix's table (`args.rs`): it has no
        // gapped block (the report prints the ungapped one only).
        gapped_karlin: if options.gapped {
            lookup_protein_params(&scoring)
        } else {
            Default::default()
        },
        gumbel: options.has_gumbel_block().then(|| {
            lookup_protein_gumbel_params(&scoring, total)
                .expect("the gap costs are in the matrix's table (`args.rs`)")
        }),
        gapped: options.gapped,
        show_linked_set_size: !options.gapped && options.sum_stats,
        // BLASTX's sections read the counts of `BlastxPairwiseOptions` (`write_blastx_pairwise_report`).
        num_descriptions: usize::MAX,
        num_alignments: usize::MAX,
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
/// The prolog of a format (`options`), whose warnings go to `diagnostics` and the posts of its
/// `-html` references where `query_warnings::html_posts_sink` sends them: `html_posts`, else
/// `diagnostics`.
pub(crate) fn write_prolog(
    writer: &mut impl crate::report::ostream::ReportWrite,
    diagnostics: &mut dyn Write,
    html_posts: Option<&mut dyn Write>,
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
) -> Result<()> {
    if options.pairwise_view().is_some() {
        crate::report::pairwise::write_blastx_prolog(
            writer,
            options.program().upper,
            &pairwise_metadata(options, subjects, shown),
            shown.pairwise_database(),
            options.html.as_deref().map(|links| {
                let posts = crate::report::query_warnings::html_posts_sink(
                    html_posts.map(|posts| posts as &mut dyn Write),
                    &mut *diagnostics,
                );
                (links, posts)
            }),
        )?;
    }
    // Outfmt 20's line of the field names, after its formatter's warning
    // (`tabular::write_prolog`).
    if let Some(spec) = options.tabular {
        super::tabular::write_prolog(
            options.program(),
            spec,
            options.taxdb_found,
            &mut *writer,
            diagnostics,
        )?;
    }
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2223-2245
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2264-2288
pub(crate) fn write_epilog(
    writer: &mut impl crate::report::ostream::ReportWrite,
    diagnostics: &mut dyn Write,
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    queries: usize,
    xml2: Option<&mut super::blast_xml::RunXml2>,
    xml_matrix: Option<crate::blastinput::query_batch::XmlFailingQuery>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2204-2212 (cited at
    // `blast_xml::matrix::XmlMatrixFailure`)
    // A run of one query whose XML2 or JSON report is held makes its report data in
    // `PrintEpilog`: a matrix name that NCBI's XML does not know fails there, before the epilog
    // (`xml_matrix`, the run's `blast_xml::xml_matrix_failure`; the XML formats of the run).
    if matches!(options.outfmt, 5 | 13..=16)
        && xml_matrix.is_some_and(|failure| failure.fails_at_end(queries))
    {
        return Err(crate::report::blast_xml::matrix::unsupported_score_matrix());
    }
    if options.pairwise_view().is_some() {
        let mut report = pairwise_metadata(options, subjects, shown);
        // NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:191-223
        // RPSTBLASTN's options name a matrix once a query batch's search has set it up
        // (`CreateRpsStructures`; every batch holds a query): a run without a query batch has a
        // null matrix name (`pairwise::write_matrix_line`).
        if options.rpstblastn.is_some() && queries == 0 {
            report.matrix_name = None;
        }
        crate::report::pairwise::write_blastx_epilog(
            writer,
            &report,
            &crate::report::pairwise::BlastxPairwiseOptions {
                word_threshold: options.threshold,
                gapped: options.gapped,
                sum_stats: options.sum_stats,
                num_descriptions: options.num_descriptions as usize,
                num_alignments: options.num_alignments as usize,
            },
            shown.pairwise_database(),
            options.html.is_some(),
        )?;
    } else if options.outfmt == 5 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2242-2243
        super::blast_xml::write_xml_epilog(
            writer,
            diagnostics,
            options.program().lower,
            queries,
            None,
        )?;
    } else if xml2.is_some() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2203-2219 and
        // :2419-2423, :2440-2445: the report a `-subject` run holds back, then the end of the
        // document.
        super::blast_xml::write_xml_epilog(
            writer,
            diagnostics,
            options.program().lower,
            queries,
            xml2,
        )?;
    } else if options.outfmt == 7 {
        // The warning of the epilog's formatter, for a field list that shows a plural name of
        // the taxonomy database without `taxdb.bti` (`tabular::write_epilog_warning`).
        if let Some(spec) = options.tabular {
            super::tabular::write_epilog_warning(
                options.program(),
                spec,
                options.taxdb_found,
                &mut *writer,
                diagnostics,
            )?;
        }
        writeln!(writer, "# BLAST processed {queries} queries")?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
        writer.ncbi_flush()?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:534-543
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:282-292
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:37-40
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2815-2822
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:94-101
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:634-639
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1042-1052
// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:294-303
// Global errorId0 "Could..." precedes per-query "Sequence..." after Combine.
pub(super) fn query_warning(
    options: &ResolvedOptions,
    query: &FastaRecord,
    search_skipped: bool,
) -> Vec<u8> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:173
    // RPSTBLASTN's results of a database of several volumes are merged, which keeps the
    // messages of error severity only (`ResolvedOptions::merged_volumes`).
    if options.merged_volumes() {
        return Vec::new();
    }
    let empty_database = options.filtered_empty && !search_skipped;
    if !query.sequence.is_empty() && !search_skipped && !empty_database {
        return Vec::new();
    }
    // The id and the title, cut as every program's query warning (`RecordNames::warning_query_id`).
    let id = query.names().warning_query_id();
    let mut bytes = format!("Warning: [{}] ", options.program().lower).into_bytes();
    bytes.extend_from_slice(&id);
    bytes.extend_from_slice(b": ");
    if search_skipped {
        bytes.extend_from_slice(b"Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options ");
    }
    if query.sequence.is_empty() {
        bytes.extend_from_slice(b"Sequence contains no data ");
    }
    if empty_database {
        bytes.extend_from_slice(EMPTY_DATABASE_MESSAGE.as_bytes());
        bytes.push(b' ');
    }
    bytes.push(b'\n');
    bytes
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:123-135
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1036-1042
#[cfg(test)]
#[path = "report_tests.rs"]
mod tests;
