//! BLASTX's tabular reports with a field list or a delimiter (outfmt 6 and 7 other than the
//! default `std` fields with a tab) and outfmt 10 and 20: the adapter of the shared writer
//! (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query of a query batch in order it
//! writes the query's warnings, outfmt 7's comment lines and one row per HSP of the query's final
//! hit list that has a Seq-align (`report::has_seqalign`), in the list's order; outfmt 20's
//! header line comes with the prolog (`write_prolog`, before the search) and outfmt 7's count of
//! queries with the epilog (`report::write_epilog`).
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
use super::input::FastaRecord;
use super::report::{
    alignment, has_seqalign, local_id, query_endpoints, query_warning, subject_label,
    tabular_io_error, ShownSubjects, TabularFlushError,
};
use super::runtime::BatchResults;
use crate::report::tabular::{
    coverage_percent, hsp_coverage, taxdb_warning_line, union_length, DeflineIds, TabularField,
    TabularHsp, TabularKind, TabularQuery, TabularSpec, TabularSubject, TabularTaxonomy,
    TabularWriter,
};

/// outfmt 7's `# BLASTX 2.17.0+` line.
const PROGRAM_VERSION: &str = concat!("BLASTX ", crate::ncbi_blast_version!());

/// The program of NCBI's messages.
const PROGRAM: &str = "blastx";

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
// ```c++
// void CBlastTabularInfo::x_CheckTaxDB()
// {
//    	if( x_IsFieldRequested(eSubjectSciNames) ||
//     	x_IsFieldRequested(eSubjectCommonNames) ||
//     	x_IsFieldRequested(eSubjectBlastNames) ||
//     	x_IsFieldRequested(eSubjectSuperKingdoms))
//    	{
//    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
//    		if(resolved.empty())
//    			ERR_POST(Warning << "Taxonomy name lookup from taxid requires installation of taxdb database with ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz");
//    	}
// }
// ```
/// The warning that each formatter of `spec` writes when it is made (`x_CheckTaxDB`), after
/// the report's lines before it are flushed; nothing when `taxdb.bti` is found or no plural name
/// is among the fields.
fn write_taxdb_warning(
    spec: &TabularSpec,
    taxdb_found: bool,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.checks_taxdb() && !taxdb_found {
        writer.flush()?;
        diagnostics.write_all(taxdb_warning_line(PROGRAM).as_bytes())?;
        diagnostics.flush()?;
    }
    Ok(())
}

/// Writes NCBI's `PrintProlog` of a custom tabular report, before the search: for outfmt 20 the
/// warning of the formatter it makes and the line of the field names; nothing for the other
/// kinds.
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
/// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254
/// ```c++
///         formatter.PrintProlog();
/// ```
/// BLASTX writes its prologs before the search, as NCBI (outfmt 0's too): a refusal of the first
/// query batch comes after them.
pub(super) fn write_prolog(
    spec: &TabularSpec,
    taxdb_found: bool,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.kind != TabularKind::CsvWithHeader {
        return Ok(());
    }
    write_taxdb_warning(spec, taxdb_found, writer, diagnostics)?;
    TabularWriter {
        spec,
        program_version: PROGRAM_VERSION,
        matrix: None,
    }
    .prolog(writer)?;
    Ok(())
}

/// The warning of the formatter that outfmt 7's `PrintEpilog` makes for its last line
/// (`report::write_epilog` writes the line).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
/// ```c++
///     if (m_FormatType == CFormattingArgs::eTabularWithComments) {
///         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec);
///         tabinfo.PrintNumProcessed(m_QueriesFormatted);
///         return;
/// ```
pub(super) fn write_epilog_warning(
    spec: &TabularSpec,
    taxdb_found: bool,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.kind == TabularKind::Commented {
        write_taxdb_warning(spec, taxdb_found, writer, diagnostics)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908
// ```c++
//     if (x_IsFieldRequested(eQueryStart) || x_IsFieldRequested(eQueryEnd) ||
//         x_IsFieldRequested(eSubjectStart) || x_IsFieldRequested(eSubjectEnd) ||
//         x_IsFieldRequested(eAlignmentLength) || x_IsFieldRequested(eGaps) ||
//         x_IsFieldRequested(eGapOpenings) || x_IsFieldRequested(eQuerySeq) ||
//         x_IsFieldRequested(eSubjectSeq) || (x_IsFieldRequested(eNumIdentical) && num_ident > 0 ) ||
//         x_IsFieldRequested(ePositives) || x_IsFieldRequested(eMismatches) ||
//         x_IsFieldRequested(ePercentPositives) || x_IsFieldRequested(ePercentIdentical) ||
//         x_IsFieldRequested(eQueryFrame) || x_IsFieldRequested(eBTOP) ||
//         x_IsFieldRequested(eSubjectStrand)) {
// ```
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
    id: &'a [u8],
    accession: &'a [u8],
    accession_version: &'a [u8],
    gi: i64,
    all_ids: std::borrow::Cow<'a, [DeflineIds]>,
    titles: Option<Vec<&'a [u8]>>,
    length: usize,
    taxonomy: &'a TabularTaxonomy,
}

/// The values of the subject `subject` (rank or index `oid`) of a report of `shown`; `slen` is the
/// length of the record, also with `-subject_loc`.
fn subject_values<'a>(
    shown: &'a ShownSubjects<'_>,
    subject: &'a FastaRecord,
    oid: usize,
    no_taxonomy: &'a TabularTaxonomy,
) -> SubjectValues<'a> {
    match shown {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:577-582
        // ```c++
        //         list<CRef<objects::CSeq_id> > subject_id_list;
        //         ITERATE(CBioseq_Handle::TId, itr, bh.GetId()) {
        //             CRef<CSeq_id> next_id = s_ReplaceLocalId(bh, itr->GetSeqId(), m_ParseLocalIds);
        //             subject_id_list.push_back(next_id);
        //         }
        //         m_SubjectIds.push_back(subject_id_list);
        // ```
        // A FASTA subject has no definition lines: its `sall*` are its own ids (the first
        // word of its title), its title and taxonomy are empty (`N/A`).
        ShownSubjects::Fasta { parsed, .. }
            if parsed.tabular.get(oid).is_some_and(Option::is_some) =>
        {
            let ids = parsed.tabular[oid].as_ref().expect("the guard");
            SubjectValues {
                id: ids.ids[0].as_bytes(),
                accession: ids.ids[1].as_bytes(),
                accession_version: ids.ids[2].as_bytes(),
                gi: parsed.subjects[oid]
                    .as_deref()
                    .map_or(0, |subject| subject.tabular_gi()),
                all_ids: std::borrow::Cow::Borrowed(&ids.all_ids),
                titles: None,
                length: subject.placement().length,
                taxonomy: no_taxonomy,
            }
        }
        // A record without a title: its `sseqid` is the first word of its generated title
        // (`report::subject_label`), its `sall*` its local id (`s_ReplaceLocalId`).
        ShownSubjects::Fasta { .. } => {
            let id = subject_label(subject).0.as_bytes();
            SubjectValues {
                id,
                accession: id,
                accession_version: id,
                gi: 0,
                all_ids: if subject.title.is_empty() {
                    std::borrow::Cow::Owned(vec![DeflineIds {
                        id: subject.internal_id.as_bytes().to_vec(),
                        gi: 0,
                        accession: subject.internal_id.as_bytes().to_vec(),
                    }])
                } else {
                    std::borrow::Cow::Borrowed(&[])
                },
                titles: None,
                length: subject.placement().length,
                taxonomy: no_taxonomy,
            }
        }
        // A database sequence: the default report's ids (`ReportedSubjects::tabular_ids`), with
        // parsed Seq-ids the ids and titles of each definition line, else its one stored title.
        ShownSubjects::Database { reported, .. } => {
            let [id, accession, accession_version] = reported.tabular_ids(oid);
            let parsed = reported.parsed_subject(oid);
            SubjectValues {
                id: id.as_bytes(),
                accession: accession.as_bytes(),
                accession_version: accession_version.as_bytes(),
                gi: parsed.map_or(0, |subject| subject.tabular_gi()),
                all_ids: std::borrow::Cow::Borrowed(reported.tabular_all_ids(oid)),
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
/// stream failure of NCBI's formatter (`report::TabularFlushError`).
#[allow(clippy::too_many_arguments)]
pub(super) fn write_blastx_tabular(
    writer: &mut impl Write,
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
    let needs = spec.needs();
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:183
    // ```c++
    //     CAlignFormatUtil::GetAsciiProteinMatrix(m_MatrixName, m_ScoringMatrix);
    // ```
    // The search's matrix (`-matrix`, item SC6 of line NP).
    let table = TabularWriter {
        spec,
        program_version: PROGRAM_VERSION,
        matrix: Some(options.scoring),
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
    // ```c++
    //         alnVec->SetGapChar('-');
    //         alnVec->SetGenCode(m_QueryGeneticCode, 0);
    //         alnVec->SetGenCode(m_DbGeneticCode, 1);
    //         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
    //         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
    // ```
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1305-1316
    // ```c++
    //     m_Ostream << "# Database: " << dbname << "\n";
    // ```
    let database: Vec<u8> = match shown {
        ShownSubjects::Fasta { path, .. } => {
            format!("User specified sequence set (Input: {})", path.display()).into_bytes()
        }
        ShownSubjects::Database { name, .. } => name.as_bytes().to_vec(),
    };
    let names = match shown {
        ShownSubjects::Database { reported, .. } => reported.tabular_names(),
        ShownSubjects::Fasta { .. } => None,
    };
    let counts_by_fields = COUNTS_FIELDS.iter().any(|field| spec.has(*field));
    let no_taxonomy = TabularTaxonomy::default();
    let mut query_ranges: Vec<(u64, u64)> = Vec::new();
    let mut next_hsp = first_hsp;
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
            // ```c++
            //     if (results.HasWarnings()) {
            //         ERR_POST(Warning << results.GetWarningStrings());
            //     }
            // ```
            stderr.write_all(&query_warning(
                query,
                batch.search_skipped,
                options.filtered_empty,
            ))?;
            stderr.write_all(shown.query_fetches(batch.query_ordinal + qidx))?;
            // NCBI makes a formatter for each query's report, and it warns when the field list
            // asks for a plural name of the taxonomy database that it cannot find
            // (tabular.cpp:139-158, `x_CheckTaxDB`).
            write_taxdb_warning(spec, options.taxdb_found, &mut *writer, &mut *stderr)?;
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
            // ```c++
            // CBlastTabularInfo::~CBlastTabularInfo()
            // {
            //     m_Ostream.flush();
            // }
            // ```
            // As the default report's: a stream failure unwinds through the query's formatter.
            (|| -> Result<()> {
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1277-1278
                // ```c++
                //     if (align_set) {
                //        int num_hits = align_set->Get().size();
                // ```
                // A query of a batch that NCBI did not search has no Seq-align set: no count.
                let count: usize = lists
                    .iter()
                    .map(|list| {
                        list.hsps
                            .iter()
                            .filter(|hsp| has_seqalign(hsp, options.gapped))
                            .count()
                    })
                    .sum();
                let mut buf: Vec<u8> = Vec::new();
                table.query_header(
                    &mut buf,
                    query.query_line().as_bytes(),
                    &database,
                    (!batch.search_skipped).then_some(count),
                )?;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:787-790
                // ```c++
                //             if(m_QueryRange.NotEmpty())
                //             	m_QueryLength = m_QueryRange.GetLength();
                //             else
                //             	m_QueryLength = query_bh.GetBioseqLength();
                // ```
                // The searched letters: the record, or its `-query_loc` interval (`qcovhsp` below);
                // `m_QueryRange` is the range as given (`query_length`).
                let searched_length = query.sequence.len();
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:880-887
                // ```c++
                //     if (x_IsFieldRequested(eQueryLength) && m_QueryLength == 0) {
                //         ...
                //     	if(m_QueryRange.NotEmpty())
                //     		m_QueryLength = m_QueryRange.GetLength();
                //     	else
                //     		m_QueryLength = sequence::GetLength(align.GetSeq_id(0), &scope);
                // ```
                // `qlen` and the `qcovs` of the subjects: the length of the `-query_loc` range as
                // given (even where it ends past the record), else the query's.
                let query_length = options
                    .query_range
                    .map_or(searched_length, |range| range.length());
                let tabular_query = TabularQuery::of_ids(
                    query.believed(),
                    local_id(query).as_bytes(),
                    query_length,
                );
                for list in lists {
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
                    let values = subject_values(shown, subject, list.oid, &no_taxonomy);
                    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:629-640
                    // ```c++
                    //        	if((check_type & eQueryCovPerSubj)) {
                    //        		int master_coverage = align_format::CAlignFormatUtil::GetMasterCoverage(tmp_align_set);
                    //        		if (master_coverage)
                    //         	{
                    //            		double subj_coverage = 100.0 * (double) master_coverage/ (double) query_len;
                    // ```
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2969-2981
                    // ```c++
                    // int CAlignFormatUtil::GetMasterCoverage(const CSeq_align_set& alnset)
                    // {
                    //     ...
                    //         CRange<TSeqPos> seq_range = (*iter)->GetSeqRange(0);
                    //         //for minus strand
                    //         if(seq_range.GetFrom() > seq_range.GetTo()){
                    //             seq_range.Set(seq_range.GetTo(), seq_range.GetFrom());
                    // ```
                    // The HSPs of one subject are one group of `InsertSubjectScores` (a hit list
                    // holds a subject once): the union of their query intervals in the record's
                    // nucleotides, in percent of the query's length.
                    let qcovs = if needs.qcovs {
                        query_ranges.clear();
                        query_ranges.extend(shown_hsps().map(|hsp| {
                            let raw = &hsp.hsp;
                            let (start, end) =
                                query_endpoints(searched_length, raw.frame, raw.q_start, raw.q_end);
                            (start.min(end) as u64 - 1, start.max(end) as u64 - 1)
                        }));
                        coverage_percent(union_length(&mut query_ranges), query_length)
                    } else {
                        -1
                    };
                    let tabular_subject = TabularSubject {
                        id: values.id,
                        accession: values.accession,
                        accession_version: values.accession_version,
                        gi: values.gi,
                        all_ids: &values.all_ids,
                        length: values.length,
                        titles: values.titles.as_deref(),
                        taxonomy: values.taxonomy,
                        names,
                        nucleotide: false,
                        qcovs,
                        qcovus: -1,
                    };
                    for hsp in shown_hsps() {
                        // The texts of the row's lookups of the taxonomy database
                        // (`x_SetTaxInfoAll`, `x_SetTaxInfo`) follow the rows before it.
                        if let Some(names) = names {
                            let text = names.row_messages(spec, values.taxonomy);
                            if !text.is_empty() {
                                writer.write_all(&buf)?;
                                buf.clear();
                                writer.flush()?;
                                stderr.write_all(&text)?;
                                stderr.flush()?;
                            }
                        }
                        let a = alignment(batch, hsp, subject, query, &code, options.scoring)?;
                        let raw = &hsp.hsp;
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
                            evalue: hsp.evalue,
                            bit_score: hsp.bit_score,
                            score: raw.score,
                            length: a.qseq.len(),
                            gap_opens: a.gap_opens,
                            gaps: a.gaps,
                            num_ident: hsp.identity,
                            q_frame: i32::from(a.frame),
                            s_frame: 0,
                            rows,
                            qcovhsp: query_coverage_per_hsp(
                                options.gapped,
                                shown_count,
                                raw.q_end - raw.q_start,
                                (a.qstart - a.qend).abs() + 1,
                                searched_length,
                                query_length,
                            ),
                        };
                        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1102-1108
                        // ```c++
                        // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
                        //     // Add tab in front of field, except for the first field.
                        //     if (iter != m_FieldsToShow.begin())
                        //         m_Ostream << m_FieldDelimiter;
                        //     x_PrintField(*iter);
                        // }
                        // m_Ostream << "\n";
                        // ```
                        // One printed row is one HSP; the probe marks it without changing it.
                        match probe.as_deref_mut() {
                            Some(probe) => {
                                writer.write_all(&buf)?;
                                buf.clear();
                                writer.flush()?;
                                probe.begin(next_hsp);
                                table.row(&mut buf, &tabular_query, &tabular_subject, &row);
                                writer.write_all(&buf)?;
                                buf.clear();
                                writer.flush()?;
                                probe.end(next_hsp);
                            }
                            None => table.row(&mut buf, &tabular_query, &tabular_subject, &row),
                        }
                        next_hsp += 1;
                    }
                }
                writer.write_all(&buf)?;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
                // ```c++
                // CBlastTabularInfo::~CBlastTabularInfo()
                // {
                //     m_Ostream.flush();
                // }
                // ```
                writer.flush().map_err(TabularFlushError)?;
                Ok(())
            })()
            .map_err(tabular_io_error)?;
        }
    }
    Ok(())
}

/// `qcovhsp` of an HSP whose query span is `span` residues of its frame and `nucleotides` letters
/// of the record, in a subject with `hsps` HSPs that have a Seq-align. `searched_length` is the
/// length of the searched letters (the `-query_loc` interval), `query_length` that of the tabular
/// `qlen` (the range as given).
///
/// A gapped search: the HSP's own Seq-align carries `hsp_percent_coverage`, which
/// `BLASTHspListToSeqAlign` computes for a translated query against a third of the query's
/// letters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1497-1499
/// ```c++
///         if  (Blast_QueryIsTranslated(program)) {
///         	s_AddScoresToSeqAlign(seqalign, hsp, emptyList, query_length/3);
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1034-1043
/// ```c
/// double Blast_HSPGetQueryCoverage(const BlastHSP* hsp, Int4 query_length)
/// {
/// 	double pct = 0;
///     if(query_length > 0) {
///             pct = 100.0 * (double) (hsp->query.end - hsp->query.offset)/ (double) query_length;
///             if(pct < 99)
///             	pct +=0.5;
///     }
///     return pct;
/// }
/// ```
/// An ungapped search makes one Seq-align of a subject's HSPs, each a Std-seg with its scores,
/// against the query's letters (`x_UngappedHSPToStdSeg`). The formatter splits a Seq-align of
/// several into one per HSP that takes the Std-seg's scores; one of a single HSP stays as it is,
/// its scores on the Std-seg, so that the Seq-align has no `hsp_percent_coverage` and
/// `x_SetQueryCovSeqalign` counts the Seq-align's span of the query (its nucleotides) against the
/// query's length, which is the range as given with `-query_loc` (tabular.cpp:786-791).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1381-1384
/// ```c++
///     CSeq_align::TScore& score_list = retval->SetScores();
///     s_BuildScoreList(hsp, score_list, seqid_list, query_length);
///
///     return retval;
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3162-3185
/// ```c++
///         if(seg.Which() == CSeq_align::C_Segs::e_Std){
///             if(seg.GetStd().size() > 1){
///                 //has more than one stdseg. Need to seperate as each
///                 //is a distinct HSP
///                 ITERATE (CSeq_align::C_Segs::TStd, iterStdseg, seg.GetStd()){
///                     CRef<CSeq_align> aln(new CSeq_align);
///                     if((*iterStdseg)->IsSetScores()){
///                         aln->SetScore() = (*iterStdseg)->GetScores();
///                     }
///                     aln->SetSegs().SetStd().push_back(*iterStdseg);
///                     alnSetRef->Set().push_back(aln);
///                 }
///
///             } else {
///                 alnSetRef->Set().push_back(*iter);
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:737-747
/// ```c++
/// void CBlastTabularInfo::x_SetQueryCovSeqalign(const CSeq_align & align, int query_len)
/// {
/// 	double tmp = 0;
/// 	if(!align.GetNamedScore("hsp_percent_coverage", tmp)) {
/// 		int len = abs((int) (align.GetSeqStop(0) - align.GetSeqStart(0))) + 1;
/// 		tmp = 100.0 * len/(double) query_len;
/// 		if(tmp  < 99)
/// 			tmp +=0.5;
/// 	}
/// 	m_QueryCovSeqalign = (int)tmp;
/// }
/// ```
fn query_coverage_per_hsp(
    gapped: bool,
    hsps: usize,
    span: i32,
    nucleotides: i32,
    searched_length: usize,
    query_length: usize,
) -> i32 {
    let searched_length = searched_length as i32;
    if gapped {
        hsp_coverage(span, searched_length / 3)
    } else if hsps > 1 {
        hsp_coverage(span, searched_length)
    } else {
        hsp_coverage(nucleotides, query_length as i32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI blast_seqalign.cpp:1497-1499 (a third of the query's letters), :1381-1384 and
    // showalign.cpp:3162-3185 with tabular.cpp:737-747 (ungapped: the query's letters, or the
    // span's nucleotides for a subject's single HSP). Values of tab.bx.cov.6 (qn2: 135 letters).
    #[test]
    fn hsp_coverage_follows_ncbis_seq_aligns() {
        assert_eq!(query_coverage_per_hsp(true, 2, 6, 18, 135, 135), 13);
        assert_eq!(query_coverage_per_hsp(true, 1, 15, 45, 135, 135), 33);
        assert_eq!(query_coverage_per_hsp(false, 2, 6, 18, 135, 135), 4);
        assert_eq!(query_coverage_per_hsp(false, 1, 6, 18, 135, 135), 13);
        assert_eq!(query_coverage_per_hsp(true, 1, 140, 420, 420, 420), 100);
        // With `-query_loc`: the interval for the engine's coverage, the range as given for the
        // formatter's.
        assert_eq!(query_coverage_per_hsp(true, 1, 15, 45, 135, 300), 33);
        assert_eq!(query_coverage_per_hsp(false, 2, 6, 18, 135, 300), 4);
        assert_eq!(query_coverage_per_hsp(false, 1, 6, 18, 135, 300), 6);
    }
}
