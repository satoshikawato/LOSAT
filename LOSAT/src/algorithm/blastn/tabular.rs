//! BLASTN's tabular reports with a field list or a delimiter (outfmt 6 and 7 other than the
//! default `std` fields with a tab) and outfmt 10 and 20: the adapter of the shared writer
//! (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query in order it writes the query's
//! warnings, outfmt 7's comment lines and one row per HSP of the query's final hit list, in the
//! list's order; outfmt 20's header line comes before the first query of the run and outfmt 7's
//! count of queries after the last.
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
//!
//!             	if(scores != CBlastFormatUtil::eNoQuerySubjCov)
//!             		CBlastFormatUtil::InsertSubjectScores (copy_aln_set, bhandle, m_QueryRange, (CBlastFormatUtil::ESubjectScores) scores);
//!            	}
//!             tabinfo.SetQueryGeneticCode(m_QueryGenCode);
//!             tabinfo.SetDbGeneticCode(m_DbGenCode);
//!             ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
//!                     const CSeq_align& s = **itr;
//!                     tabinfo.SetFields(s, *m_Scope, &m_ScoringMatrix);
//!                     tabinfo.Print();
//!             }
//! ```
//! What `SetFields` reads of an HSP is its Seq-align (`algorithm::blastn::blast_xml` has the
//! coordinates): the query on its plus strand, the subject descending for a minus-strand HSP
//! (frames `1/-1`), the score, e-value and bit score, the alignment's lengths and gaps, and for
//! `qseq`, `sseq`, `btop`, `positive` and `ppos` the two rows, from which it counts the
//! identities again (`report::tabular::TabularWriter::row`). BLASTN's formatter reads no other
//! sequence data (`SetNoFetch(true)`) and has no scoring matrix: positives are identities.

use std::borrow::Cow;
use std::io::Write;
use std::sync::Arc;

use anyhow::{Context, Result};

use super::hsp::{blastn_hsp_count, BlastnHitList, BlastnHsp, NCBI_BLASTN_VERSION};
use super::pairwise::xml_rows_into;
use super::x_subject_store::SubjectLetters;
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::common::GapEditOp;
use crate::report::db_subject::ParsedSubject;
use crate::report::query_warnings::QueryWarnings;
use crate::report::tabular::{
    coverage_percent, hsp_coverage, taxdb_warning_line, union_length, DeflineIds, TabularField,
    TabularHsp, TabularKind, TabularNames, TabularQuery, TabularSpec, TabularSubject,
    TabularTaxonomy, TabularWriter,
};

const PROGRAM: &str = "BLASTN";

/// What the report needs besides the final hit lists.
pub(crate) struct BlastnTabular<'a> {
    pub spec: &'a TabularSpec,
    /// The queries of this call (their searched letters, which the rows read).
    pub queries: &'a [FastaRecord],
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
    /// The letters of the subjects, which the rows read.
    pub letters: SubjectLetters<'a>,
    /// By query: whether NCBI did not search its batch (outfmt 7's lines without a count).
    pub unsearched: &'a [bool],
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
pub(crate) enum TabularSubjects<'a> {
    /// The records of `-subject`: the tabular id of each (`sseqid`, `sacc`, `saccver`; the
    /// default report's) and the length of its record.
    Records {
        ids: &'a [Arc<[u8]>],
        length: &'a dyn Fn(usize) -> usize,
    },
    /// The reported sequences of a database: their tabular ids, ordinal ids, records (with
    /// the stored title of a sequence without parsed Seq-ids), parsed Seq-ids (empty when none
    /// has them), lengths, and taxonomy (empty when no field shows it).
    Database {
        ids: &'a [Arc<[u8]>],
        ordinals: &'a [u32],
        records: &'a [FastaRecord],
        parsed: &'a [Option<Arc<ParsedSubject>>],
        lengths: &'a [usize],
        taxonomy: &'a [TabularTaxonomy],
    },
}

/// The values of one subject that its rows show (`TabularSubject` without the coverages).
struct SubjectValues<'a> {
    id: Cow<'a, [u8]>,
    accession: Cow<'a, [u8]>,
    accession_version: Cow<'a, [u8]>,
    gi: i64,
    all_ids: Vec<DeflineIds>,
    titles: Option<Vec<&'a [u8]>>,
    length: usize,
    taxonomy: Cow<'a, TabularTaxonomy>,
}

impl<'a> TabularSubjects<'a> {
    /// The values of subject `s_idx` that `needs` asks for.
    fn values(&self, s_idx: usize, spec: &TabularSpec) -> Result<SubjectValues<'a>> {
        let needs = spec.needs();
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
                let id: &[u8] = ids.get(s_idx).map_or(&b"unknown"[..], |id| id.as_ref());
                Ok(SubjectValues {
                    id: Cow::Borrowed(id),
                    accession: Cow::Borrowed(id),
                    accession_version: Cow::Borrowed(id),
                    gi: 0,
                    all_ids: Vec::new(),
                    titles: None,
                    length: if spec.has(TabularField::SubjectLength) {
                        length(s_idx)
                    } else {
                        0
                    },
                    taxonomy: Cow::Owned(TabularTaxonomy::default()),
                })
            }
            TabularSubjects::Database {
                ids,
                ordinals,
                records,
                parsed,
                lengths,
                taxonomy,
            } => {
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
                                ordinals.get(s_idx).copied().unwrap_or_default(),
                                &format!("has {problem}"),
                                PROGRAM,
                            )
                        };
                        let [id, accession, accession_version] =
                            subject.tabular_ids(false).map_err(refused)?;
                        let all_ids = if needs.all_deflines {
                            subject.all_tabular_ids(false).map_err(refused)?
                        } else {
                            Vec::new()
                        };
                        Ok(SubjectValues {
                            id: Cow::Owned(id),
                            accession: Cow::Owned(accession),
                            accession_version: Cow::Owned(accession_version),
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
                    // A sequence without parsed Seq-ids: the first word of its stored title
                    // (`gnl|BL_ORD_ID|n` replaced, `CShowBlastDefline::GetSeqIdList`) for every
                    // id, and its one definition line's title.
                    None => {
                        let id: &[u8] = ids.get(s_idx).map_or(&b"unknown"[..], |id| id.as_ref());
                        let title = records
                            .get(s_idx)
                            .map_or(&[][..], |record| record.title.as_slice());
                        Ok(SubjectValues {
                            id: Cow::Borrowed(id),
                            accession: Cow::Borrowed(id),
                            accession_version: Cow::Borrowed(id),
                            gi: 0,
                            all_ids: Vec::new(),
                            titles: Some(vec![title]),
                            length,
                            taxonomy,
                        })
                    }
                }
            }
        }
    }
}

/// Writes outfmt 20's line of the field names (`PrintProlog`; nothing for the other kinds).
/// The taxonomy database's warning of the formatter that NCBI makes for it is written before the
/// search (`write_tabular_prolog_warning`).
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
pub(crate) fn write_tabular_prolog(spec: &TabularSpec, writer: &mut dyn Write) -> Result<()> {
    let version = format!("BLASTN {NCBI_BLASTN_VERSION}");
    TabularWriter {
        spec,
        program_version: &version,
        matrix: None,
    }
    .prolog(writer)?;
    Ok(())
}

/// Writes the warning of the formatter that NCBI makes for outfmt 20's header line, which it
/// writes before it reads the first query batch (`PrintProlog`, blastn_app.cpp:258): the
/// taxonomy database's warning when the field list checks it and `taxdb.bti` is not found
/// (`x_CheckTaxDB`, `TabularSpec::checks_taxdb`). A run that finds it writes nothing here, so
/// the line itself can come with the first report, after LOSAT's refusals.
pub(crate) fn write_tabular_prolog_warning(
    spec: &TabularSpec,
    taxdb_found: bool,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.kind == TabularKind::CsvWithHeader && spec.checks_taxdb() && !taxdb_found {
        diagnostics.write_all(taxdb_warning_line("blastn").as_bytes())?;
        diagnostics.flush()?;
    }
    Ok(())
}

/// Whether `taxdb.bti` is on NCBI's search path for a `-subject` search (`x_CheckTaxDB`'s
/// `SeqDB_ResolveDbPath("taxdb.bti")`: the current directory, `BLASTDB` and the registry's
/// `[BLAST] BLASTDB`, as a database is looked for, `seqdb::SearchPath`). A registry value that
/// LOSAT does not read as NCBI is its refusal, as for `-db`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
/// ```c++
/// void CBlastTabularInfo::x_CheckTaxDB()
/// {
///    	if( x_IsFieldRequested(eSubjectSciNames) ||
///     ...
///    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
///    		if(resolved.empty())
/// ```
pub(crate) fn subject_taxdb_found(
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<bool> {
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => {
            anyhow::bail!("{reason}; this is not supported by LOSAT's {PROGRAM}")
        }
    };
    let path = crate::seqdb::SearchPath::from_environment(registry)
        .map_err(|error| crate::algorithm::common::database::database_error(error, PROGRAM))?;
    Ok(path.resolve_file("taxdb.bti").is_some())
}

/// Writes the report of every query of `hit_lists` (indexed by query) to `writer`, each query's
/// warnings (and the taxonomy database's warning of its formatter) before its lines, outfmt 20's
/// header line first when `tab.prolog` and outfmt 7's last line after the last query when
/// `tab.epilog`. `probe` marks each row (the format observer).
pub(crate) fn write_blastn_tabular(
    tab: &BlastnTabular<'_>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut dyn Write,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    mut probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    let spec = tab.spec;
    let needs = spec.needs();
    // The refusals of `sall*` (`ParsedSubject::all_tabular_ids`) come before the first byte;
    // those of the other ids are the reader's (`ShownSubjects`).
    for hit_list in hit_lists.iter().flatten().filter(|_| needs.all_deflines) {
        for hsp_list in &hit_list.hsplist_array {
            if let Some(first) = hsp_list.hsps.first() {
                tab.subjects.values(first.s_idx as usize, spec)?;
            }
        }
    }
    let version = format!("BLASTN {NCBI_BLASTN_VERSION}");
    let table = TabularWriter {
        spec,
        program_version: &version,
        matrix: None,
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
        (spec.checks_taxdb() && !tab.taxdb_found).then(|| taxdb_warning_line("blastn"));
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
    let (mut query_row, mut subject_row): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    let mut query_ranges: Vec<(u64, u64)> = Vec::new();
    let mut hsp_index: HspIndex = 0;
    for (q_idx, hit_list) in hit_lists.iter().enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`), then the warning of its
        // formatter.
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        warn(&mut *writer, warnings.as_deref_mut())?;
        let hit_list = hit_list.as_ref();
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
            id: tab
                .query_ids
                .get(q_idx)
                .map_or(&b"unknown"[..], |id| id.as_ref()),
            length: query_length,
        };
        buf.clear();
        // A query of a batch that NCBI did not search has no count (`write_query_header`).
        table.query_header(
            &mut buf,
            tab.query_titles
                .get(q_idx)
                .map_or(&b"unknown"[..], |title| title.as_ref()),
            tab.database.as_bytes(),
            (!tab.unsearched.get(q_idx).copied().unwrap_or(false))
                .then(|| blastn_hsp_count(hit_list)),
        )?;
        if probe.is_some() {
            writer.write_all(&buf)?;
            buf.clear();
        }
        let Some(hit_list) = hit_list else {
            writer.write_all(&buf)?;
            continue;
        };
        for hsp_list in &hit_list.hsplist_array {
            let Some(first) = hsp_list.hsps.first() else {
                continue;
            };
            let s_idx = first.s_idx as usize;
            let values = tab.subjects.values(s_idx, spec)?;
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
            // subject once).
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
            let qcovus = if needs.qcovus {
                uniq_subject_coverage(&hsp_list.hsps, query_length)
            } else {
                -1
            };
            let subject = TabularSubject {
                id: &values.id,
                accession: &values.accession,
                accession_version: &values.accession_version,
                gi: values.gi,
                all_ids: &values.all_ids,
                length: values.length,
                titles: values.titles.as_deref(),
                taxonomy: &values.taxonomy,
                names: tab.names,
                nucleotide: true,
                qcovs,
                qcovus,
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
                if needs.rows {
                    xml_rows_into(
                        hsp,
                        tab.queries,
                        tab.letters,
                        &mut query_row,
                        &mut subject_row,
                    )
                    .context("BLASTN tabular rows of a final HSP")?;
                }
                let row = tabular_hsp(
                    hsp,
                    searched_length,
                    needs
                        .rows
                        .then_some((query_row.as_slice(), subject_row.as_slice())),
                );
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

/// The values of one HSP's row: its Seq-align's coordinates, scores and lengths (`SetFields`),
/// `rows` when a field shows them.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1061-1091
/// ```c++
///         if(x_IsFieldRequested(eSubjectStrand))
///         {
///         	if(!subject_is_na)
///         		m_SubjectStrand = NA;
///         	else
///         		m_SubjectStrand = ((s_start - s_end) > 0 )? "minus":"plus";
///         }
///     }
///     SetEndpoints(q_start, q_end, s_start, s_end);
///
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
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1100-1121
/// ```c++
///             if (chunk->IsGap()) {
///                 ++num_gap_opens;
///                 num_gaps += chunk_length;
///             }
/// ```
/// The gaps are the letters of the edit script's gap operations; the gap openings and the
/// identities are the default row's (`write_hsp_rows`).
fn tabular_hsp<'r>(
    hsp: &BlastnHsp,
    searched_length: usize,
    rows: Option<(&'r [u8], &'r [u8])>,
) -> TabularHsp<'r> {
    let gaps = hsp.gap_info.as_deref().map_or(0, |ops| {
        ops.iter()
            .filter(|op| !matches!(op, GapEditOp::Sub(_)))
            .map(|op| op.num() as usize)
            .sum()
    });
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
        gaps,
        num_ident: hsp.num_ident,
        q_frame: 1,
        s_frame: if hsp.s_start > hsp.s_end { -1 } else { 1 },
        rows,
        // `Blast_HSPGetQueryCoverage` of the HSP's context, which has the searched length.
        qcovhsp: hsp_coverage((hsp.q_end - hsp.q_start + 1) as i32, searched_length as i32),
    }
}

/// `qcovus` of one subject's HSPs (in the Seq-align order) for a query of `query_length`:
/// `GetUniqSeqCoverage` in percent as `InsertSubjectScores` keeps it, -1 (`N/A`) when nothing is
/// covered.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:646-657
/// ```c++
///        	if((check_type & eQueryCovPerUniqSubj)) {
///        		int uniq_coverage = align_format::CAlignFormatUtil::GetUniqSeqCoverage(tmp_align_set);
///        		if (uniq_coverage)
///        	    {
///        	    	double uniq_subj_coverage = 100.0 * (double) uniq_coverage/ (double) query_len;
///        	    	if(uniq_subj_coverage < 99)
///        	        	uniq_subj_coverage +=0.5;
///
///        	         (*left_it)->SetNamedScore ("uniq_seq_percent_coverage", (int) uniq_subj_coverage);
///        	    }
///        	}
/// ```
fn uniq_subject_coverage(hsps: &[BlastnHsp], query_length: usize) -> i32 {
    if query_length == 0 {
        return -1;
    }
    // `int uniq_coverage` of the `TSeqPos` length.
    let covered = uniq_seq_coverage(hsps) as i32;
    if covered == 0 {
        return -1;
    }
    let mut pct = 100.0 * f64::from(covered) / query_length as f64;
    if pct < 99.0 {
        pct += 0.5;
    }
    pct as i32
}

/// NCBI's `CRange<TSeqPos>`: `from` and `to_open` as `TSeqPos` (unsigned 32-bit, wrapping).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SeqRange {
    from: u32,
    to_open: u32,
}

impl SeqRange {
    /// `CRange(from, to)` (`Set`: `SetFrom(from).SetTo(to)`, `SetToOpen(to + 1)`).
    fn new(from: u32, to: u32) -> Self {
        SeqRange {
            from,
            to_open: to.wrapping_add(1),
        }
    }

    fn to(self) -> u32 {
        self.to_open.wrapping_sub(1)
    }

    fn is_empty(self) -> bool {
        self.to_open <= self.from
    }

    fn length(self) -> u32 {
        if self.to_open <= self.from {
            0
        } else {
            self.to_open - self.from
        }
    }

    /// `CombineWith`.
    fn combine_with(&mut self, r: SeqRange) {
        if !r.is_empty() {
            if !self.is_empty() {
                self.from = self.from.min(r.from);
                self.to_open = self.to_open.max(r.to_open);
            } else {
                *self = r;
            }
        }
    }

    /// `s_FixMinusStrandRange`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4829-4836
    /// ```c++
    /// static inline CRange<TSeqPos> & s_FixMinusStrandRange(CRange<TSeqPos> & rng)
    /// {
    /// 	if(rng.GetFrom() > rng.GetTo()){
    /// 		rng.Set(rng.GetTo(), rng.GetFrom());
    /// 	}
    /// ```
    fn fix_minus_strand(self) -> Self {
        if self.from > self.to() {
            SeqRange::new(self.to(), self.from)
        } else {
            self
        }
    }
}

/// NCBI's `CRangeCollection<TSeqPos>`: sorted ranges that do not overlap or touch (the parts of
/// `x_CombineWith`, `x_Subtract` and `GetCoveredLength` that `GetUniqSeqCoverage` calls, with
/// their `TSeqPos` arithmetic).
#[derive(Clone, Debug, Default)]
struct RangeCollection(Vec<SeqRange>);

impl RangeCollection {
    /// `std::lower_bound` from `start` with `PRangeLessPos` (`R.GetToOpen() <= Pos`): the first
    /// range with `to_open > pos` while the ranges are sorted. libstdc++'s halving, step for
    /// step, since a range that a query gap made (`to_open` 0, `uniq_seq_coverage`) leaves the
    /// collection unsorted, where the probes decide the result.
    ///
    /// libstdc++ reference (bits/stl_algobase.h, `std::__lower_bound`):
    /// ```c++
    ///       while (__len > 0)
    /// 	{
    /// 	  _DistanceType __half = __len >> 1;
    /// 	  _ForwardIterator __middle = __first;
    /// 	  std::advance(__middle, __half);
    /// 	  if (__comp(__middle, __val))
    /// 	    {
    /// 	      __first = __middle;
    /// 	      ++__first;
    /// 	      __len = __len - __half - 1;
    /// 	    }
    /// 	  else
    /// 	    __len = __half;
    /// 	}
    ///       return __first;
    /// ```
    fn lower_bound(&self, start: usize, pos: u32) -> usize {
        let mut first = start;
        let mut len = self.0.len() - start;
        while len > 0 {
            let half = len >> 1;
            let middle = first + half;
            if self.0[middle].to_open <= pos {
                first = middle + 1;
                len = len - half - 1;
            } else {
                len = half;
            }
        }
        first
    }

    /// `x_CombineWith(const TRange&)`.
    ///
    /// NCBI reference (598d8ae6): c++/include/util/range_coll.hpp:321-346
    /// ```c++
    ///         position_type pos_from = r.GetFrom();
    ///         position_type pos_to_open = r.GetToOpen();
    ///
    ///         if (pos_from == TRange::GetWholeFrom()) {
    ///             pos_from++; // Prevent overflow
    ///         }
    ///         iterator it_begin_m =
    ///             lower_bound(begin_nc(), end_nc(), pos_from - 1, p); /* NCBI_FAKE_WARNING: WorkShop */
    ///         if(it_begin_m != end_nc() && it_begin_m->GetFrom() <= pos_to_open)  { // intersection
    ///             iterator it_end_m = lower_bound(it_begin_m, end_nc(), pos_to_open, p);
    ///             it_begin_m->CombineWith(r);
    ///
    ///             if(it_end_m != end_nc()  &&  it_end_m->GetFrom() <= pos_to_open) {// subject to merge
    ///                 it_begin_m->SetToOpen(it_end_m->GetToOpen());
    ///                 ++it_end_m; // including it into erased set
    ///             }
    ///             m_vRanges.erase(it_begin_m + 1, it_end_m);
    ///         } else {
    ///             m_vRanges.insert(it_begin_m, r);
    ///         }
    /// ```
    fn combine_with(&mut self, r: SeqRange) {
        let mut pos_from = r.from;
        let pos_to_open = r.to_open;
        if pos_from == 0 {
            pos_from += 1;
        }
        let begin = self.lower_bound(0, pos_from - 1);
        if begin < self.0.len() && self.0[begin].from <= pos_to_open {
            let mut end = self.lower_bound(begin, pos_to_open);
            self.0[begin].combine_with(r);
            if end < self.0.len() && self.0[end].from <= pos_to_open {
                self.0[begin].to_open = self.0[end].to_open;
                end += 1;
            }
            self.0.drain(begin + 1..end);
        } else {
            self.0.insert(begin, r);
        }
    }

    /// `x_Subtract(const TRange&)`.
    ///
    /// NCBI reference (598d8ae6): c++/include/util/range_coll.hpp:348-376
    /// ```c++
    ///         iterator it_begin_e = lower_bound(begin_nc(), end_nc(), pos_from, P);
    ///         if(it_begin_e != end_nc()) {   // possible intersection
    ///
    ///             if(it_begin_e->GetFrom() < pos_from  &&  it_begin_e->GetToOpen() > pos_to_open)    {
    ///                 //it_begin_e contains R, split it in two
    ///                 TRange it_r = *it_begin_e;
    ///                 it_begin_e = m_vRanges.insert(it_begin_e, it_r);
    ///                 it_begin_e->SetToOpen(pos_from);
    ///                 (++it_begin_e)->SetFrom(pos_to_open);
    ///             } else  {
    ///                 if(it_begin_e->GetFrom() < pos_from) { // cut this segement , but don't erase
    ///                     it_begin_e->SetToOpen(pos_from);
    ///                     ++it_begin_e;
    ///                 }
    ///                 iterator it_end_e = lower_bound(it_begin_e, end_nc(), pos_to_open, P);
    ///                 if(it_end_e != end_nc()  &&  it_end_e->GetFrom() < pos_to_open)    {
    ///                     it_end_e->SetFrom(pos_to_open); // truncate
    ///                 }
    ///                 // erase ranges fully covered by R
    ///                 m_vRanges.erase(it_begin_e, it_end_e);
    ///             }
    ///         }
    /// ```
    fn subtract(&mut self, r: SeqRange) {
        let pos_from = r.from;
        let pos_to_open = r.to_open;
        let mut begin = self.lower_bound(0, pos_from);
        if begin == self.0.len() {
            return;
        }
        if self.0[begin].from < pos_from && self.0[begin].to_open > pos_to_open {
            let range = self.0[begin];
            self.0.insert(begin, range);
            self.0[begin].to_open = pos_from;
            self.0[begin + 1].from = pos_to_open;
        } else {
            if self.0[begin].from < pos_from {
                self.0[begin].to_open = pos_from;
                begin += 1;
            }
            let end = self.lower_bound(begin, pos_to_open);
            if end < self.0.len() && self.0[end].from < pos_to_open {
                self.0[end].from = pos_to_open;
            }
            self.0.drain(begin..end);
        }
    }

    /// `GetCoveredLength`: the sum of the lengths, in `TSeqPos`.
    fn covered_length(&self) -> u32 {
        self.0
            .iter()
            .fold(0u32, |total, range| total.wrapping_add(range.length()))
    }
}

/// The subject range of an HSP's Seq-align (`GetSeqRange(1)`, record coordinates, 0-based).
fn subject_range(hsp: &BlastnHsp) -> SeqRange {
    let (low, high) = (hsp.s_start.min(hsp.s_end), hsp.s_start.max(hsp.s_end));
    SeqRange::new(low as u32 - 1, high as u32 - 1)
}

/// The query range of an HSP's Seq-align (`GetSeqRange(0)`: the query always ascending).
fn query_range(hsp: &BlastnHsp) -> SeqRange {
    SeqRange::new(hsp.q_start as u32 - 1, hsp.q_end as u32 - 1)
}

/// The query position (record coordinates, 0-based) in the column of the HSP's alignment that
/// holds subject position `subject_pos`, -1 where the query has a gap there: `CAlnMap`'s
/// `GetSeqPosFromAlnPos(0, GetAlnPosFromSeqPos(1, subject_pos))` on the HSP's Dense-seg, whose
/// subject row ascends and whose query row descends for a minus-strand HSP
/// (`s_CollectSeqAlignData`; the columns are the edit script's).
///
/// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:663-686
/// ```c++
///     TNumseg seg = GetSeg(aln_pos);
///     TSignedSeqPos pos = GetStart(for_row, seg);
///     if (pos >= 0) {
///         TSeqPos delta = (aln_pos - GetAlnStart(seg)) * GetWidth(for_row);
///         if (IsPositiveStrand(for_row)) {
///             pos += delta;
///         } else {
///             pos += x_GetLen(for_row, x_GetRawSegFromSeg(seg)) - 1 - delta;
///         }
///     } else if (dir != eNone) {
/// ```
fn query_position(hsp: &BlastnHsp, subject_pos: u32) -> i64 {
    let subject_first = hsp.s_start.min(hsp.s_end) as i64 - 1;
    let target = i64::from(subject_pos) - subject_first;
    let ungapped = [GapEditOp::Sub(
        (hsp.internal_q_end_0 - hsp.internal_q_offset_0) as u32,
    )];
    let ops = hsp.gap_info.as_deref().unwrap_or(&ungapped);
    // The query letters and subject letters before the current operation.
    let (mut query, mut subject) = (0i64, 0i64);
    for op in ops {
        let n = i64::from(op.num());
        match op {
            GapEditOp::Sub(_) => {
                if target < subject + n {
                    let offset = query + (target - subject);
                    return if hsp.query_frame < 0 {
                        hsp.q_end as i64 - 1 - offset
                    } else {
                        hsp.q_start as i64 - 1 + offset
                    };
                }
                query += n;
                subject += n;
            }
            GapEditOp::Del(_) => {
                if target < subject + n {
                    return -1;
                }
                subject += n;
            }
            GapEditOp::Ins(_) => query += n,
        }
    }
    -1
}

/// `GetUniqSeqCoverage` of one subject's HSPs in the Seq-align order: the query letters that
/// the HSPs cover, where a part of an HSP whose subject letters an HSP before it covered does
/// not count (its query positions mapped through the alignment's columns).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4838-4892
/// ```c++
/// 	list<CRef<CSeq_align> >::iterator mItr=alnset.Set().begin();
/// 	CRangeCollection<TSeqPos> subj_rng_coll((*mItr)->GetSeqRange(1));
/// 	CRange<TSeqPos> q_rng((*mItr)->GetSeqRange(0));
/// 	CRangeCollection<TSeqPos> query_rng_coll(s_FixMinusStrandRange(q_rng));
/// 	++mItr;
/// 	for(;mItr != alnset.Set().end(); ++mItr) {
/// 		const CRange<TSeqPos> align_subj_rng((*mItr)->GetSeqRange(1));
/// 		// subject range should always be on the positive strand
/// 		ASSERT(align_subj_rng.GetTo() > align_subj_rng.GetFrom());
/// 		CRangeCollection<TSeqPos> coll(align_subj_rng);
/// 		coll.Subtract(subj_rng_coll);
///
/// 		if (coll.empty())
/// 			continue;
///
/// 		if(coll[0] == align_subj_rng) {
/// 			CRange<TSeqPos> query_rng ((*mItr)->GetSeqRange(0));
/// 			query_rng_coll += s_FixMinusStrandRange(query_rng);
/// 			subj_rng_coll += align_subj_rng;
/// 		}
/// 		else {
/// 			ITERATE (CRangeCollection<TSeqPos>, uItr, coll) {
/// 				CRange<TSeqPos> query_rng;
/// 				const CRange<TSeqPos> & subj_rng = (*uItr);
/// 				...
/// 				CAlnMap map(densegAln->GetSegs().GetDenseg());
/// 				TSignedSeqPos subj_aln_start =  map.GetAlnPosFromSeqPos(1,subj_rng.GetFrom());
/// 				TSignedSeqPos subj_aln_end =  map.GetAlnPosFromSeqPos(1,subj_rng.GetTo());
/// 				query_rng.SetFrom(map.GetSeqPosFromAlnPos(0,subj_aln_start));
/// 				query_rng.SetTo(map.GetSeqPosFromAlnPos(0,subj_aln_end));
///
/// 				query_rng_coll += s_FixMinusStrandRange(query_rng);
/// 				subj_rng_coll += subj_rng;
/// 			}
/// 		}
/// 	}
///
/// 	return query_rng_coll.GetCoveredLength();
/// ```
/// A query gap in the column of a part's first or last subject letter gives NCBI the position
/// -1, which the `TSeqPos` range keeps as 2^32 - 1 (the same arithmetic here). LOSAT's BLASTN
/// alignments are always gapped Dense-segs (no Dense-diag).
fn uniq_seq_coverage(hsps: &[BlastnHsp]) -> u32 {
    let Some((first, rest)) = hsps.split_first() else {
        return 0;
    };
    let mut subject_ranges = RangeCollection(vec![subject_range(first)]);
    let mut query_ranges = RangeCollection(vec![query_range(first).fix_minus_strand()]);
    for hsp in rest {
        let align_subject = subject_range(hsp);
        let mut parts = RangeCollection(vec![align_subject]);
        for &covered in &subject_ranges.0 {
            parts.subtract(covered);
        }
        if parts.0.is_empty() {
            continue;
        }
        if parts.0[0] == align_subject {
            query_ranges.combine_with(query_range(hsp).fix_minus_strand());
            subject_ranges.combine_with(align_subject);
        } else {
            for part in parts.0 {
                // `SetFrom` and `SetTo` of a `TSignedSeqPos` into the `TSeqPos` range.
                let from = query_position(hsp, part.from) as u32;
                let to = query_position(hsp, part.to()) as u32;
                let query_part = SeqRange {
                    from,
                    to_open: to.wrapping_add(1),
                };
                query_ranges.combine_with(query_part.fix_minus_strand());
                subject_ranges.combine_with(part);
            }
        }
    }
    query_ranges.covered_length()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collection(ranges: &[(u32, u32)]) -> RangeCollection {
        let mut out = RangeCollection::default();
        for &(from, to) in ranges {
            out.combine_with(SeqRange::new(from, to));
        }
        out
    }

    fn pairs(collection: &RangeCollection) -> Vec<(u32, u32)> {
        collection.0.iter().map(|r| (r.from, r.to())).collect()
    }

    #[test]
    fn range_collection_merges_overlapping_and_touching_ranges() {
        let merged = collection(&[(10, 19), (30, 39), (20, 25), (0, 3)]);
        assert_eq!(pairs(&merged), [(0, 3), (10, 25), (30, 39)]);
        assert_eq!(merged.covered_length(), 4 + 16 + 10);
        let spanning = collection(&[(10, 19), (30, 39), (5, 50)]);
        assert_eq!(pairs(&spanning), [(5, 50)]);
        let zero = collection(&[(5, 9), (0, 4)]);
        assert_eq!(pairs(&zero), [(0, 9)]);
    }

    // A query gap at a part's last subject letter makes the range [from, to_open 0), which is
    // empty and leaves the collection unsorted: the probes of `std::lower_bound` then decide
    // where later ranges go (`uniq_seq_coverage`).
    #[test]
    fn lower_bound_probes_as_libstdcxx() {
        let unsorted = RangeCollection(vec![
            SeqRange::new(0, 9),
            SeqRange {
                from: 1546,
                to_open: 0,
            },
            SeqRange::new(20, 29),
            SeqRange::new(40, 49),
        ]);
        // The probes: [20, 30) is below 35, then [40, 50) is not.
        assert_eq!(unsorted.lower_bound(0, 35), 3);
        // A linear scan would answer 0 ([0, 10) is not below 5); the halving probes [20, 30),
        // which is not below 5, then [1546, 0), which is, and answers 2.
        assert_eq!(unsorted.lower_bound(0, 5), 2);
        let mut collection = unsorted.clone();
        collection.combine_with(SeqRange::new(12, 14));
        assert_eq!(
            collection
                .0
                .iter()
                .map(|r| (r.from, r.to_open))
                .collect::<Vec<_>>(),
            [(0, 10), (1546, 0), (12, 15), (20, 30), (40, 50)]
        );
    }

    #[test]
    fn range_collection_subtracts_as_ncbi() {
        let mut parts = RangeCollection(vec![SeqRange::new(0, 99)]);
        parts.subtract(SeqRange::new(10, 19));
        assert_eq!(pairs(&parts), [(0, 9), (20, 99)]);
        parts.subtract(SeqRange::new(90, 120));
        assert_eq!(pairs(&parts), [(0, 9), (20, 89)]);
        parts.subtract(SeqRange::new(0, 30));
        assert_eq!(pairs(&parts), [(31, 89)]);
        parts.subtract(SeqRange::new(0, 200));
        assert!(parts.0.is_empty());
    }

    fn hsp(
        q: (usize, usize),
        s: (usize, usize),
        minus: bool,
        ops: Option<Vec<GapEditOp>>,
    ) -> BlastnHsp {
        BlastnHsp {
            identity: 0.0,
            length: 0,
            mismatch: 0,
            gapopen: 0,
            q_start: q.0,
            q_end: q.1,
            s_start: s.0,
            s_end: s.1,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: if minus { -1 } else { 1 },
            query_length: 1000,
            q_idx: 0,
            s_idx: 0,
            raw_score: 0,
            internal_q_offset_0: 0,
            internal_q_end_0: q.1 - q.0 + 1,
            internal_s_offset_0: 0,
            internal_s_end_0: s.0.max(s.1) - s.0.min(s.1) + 1,
            internal_query_context_offset: 0,
            gap_info: ops,
            num_positives: 0,
        }
    }

    #[test]
    fn query_positions_follow_the_columns() {
        // Query 11..=20 against subject 101..=111 with a query gap after 4 columns.
        let gapped = hsp(
            (11, 20),
            (101, 111),
            false,
            Some(vec![
                GapEditOp::Sub(4),
                GapEditOp::Del(1),
                GapEditOp::Sub(6),
            ]),
        );
        assert_eq!(query_position(&gapped, 100), 10);
        assert_eq!(query_position(&gapped, 103), 13);
        assert_eq!(query_position(&gapped, 104), -1);
        assert_eq!(query_position(&gapped, 105), 14);
        assert_eq!(query_position(&gapped, 110), 19);
        // A minus-strand HSP: the query descends along the subject.
        let minus = hsp((11, 20), (110, 101), true, None);
        assert_eq!(query_position(&minus, 100), 19);
        assert_eq!(query_position(&minus, 109), 10);
    }

    #[test]
    fn uniq_coverage_skips_subject_letters_covered_before() {
        // Two HSPs on disjoint subject ranges: both query ranges count.
        let a = hsp((1, 100), (1, 100), false, None);
        let b = hsp((201, 300), (201, 300), false, None);
        assert_eq!(uniq_seq_coverage(&[a.clone(), b]), 200);
        // A repeat of the query on the same subject letters counts once.
        let repeat = hsp((401, 500), (1, 100), false, None);
        assert_eq!(uniq_seq_coverage(&[a.clone(), repeat]), 100);
        // A partly covered HSP: only its uncovered subject part maps to query letters.
        let partial = hsp((501, 600), (51, 150), false, None);
        assert_eq!(uniq_seq_coverage(&[a, partial]), 100 + 50);
        assert_eq!(uniq_subject_coverage(&[], 100), -1);
    }
}
