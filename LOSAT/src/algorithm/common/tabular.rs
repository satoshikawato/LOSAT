//! The subjects of the tabular reports with a field list (outfmt 6 and 7 other than the default
//! `std` fields with a tab, outfmt 10 and 20; `report::tabular`), as the programs' adapters hand
//! them to the shared writer: the values of a FASTA subject or a database sequence that its rows
//! show, the taxonomy of the reported database sequences, and outfmt 20's line of the field names
//! with its formatter's warning. For the programs whose subjects are nucleotide sequences;
//! `program` names the program in LOSAT's refusals. BLASTN (line CLI C2, PR 1; its outfmt 20 line
//! and warning are its own, DECISIONS 6), TBLASTN and TBLASTX (PR 2) use it.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-851
//! ```c++
//!        		const CBioseq_Handle& subject_bh =
//!                 scope.GetBioseqHandle(align.GetSeq_id(1));
//!             if(setSubjectId) {
//!             	SetSubjectId(subject_bh);
//!             }
//!             m_SubjectLength = subject_bh.GetBioseqLength();
//!
//!             if(setSubjectIds || setSubjectTaxInfo || setSubjectTitle || setSubjectTaxInfoAll) {
//!             	CRef<CBlast_def_line_set> bdlRef =
//!             			CSeqDB::ExtractBlastDefline(subject_bh);
//!             	if(setSubjectIds) {
//!             		x_SetSubjectIds(subject_bh, bdlRef);
//!             	}
//! ```

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

/// The subjects as the report shows them, by subject index (a subject set) or rank (the
/// reported sequences of a database).
pub(crate) enum TabularSubjects<'a> {
    /// The records of `-subject`: the tabular id of each (`sseqid`, `sacc`, `saccver`; the
    /// default report's), the length of its record, and the ones whose ids `-parse_deflines`
    /// read (empty without the option).
    Records {
        ids: &'a [Arc<[u8]>],
        length: &'a dyn Fn(usize) -> usize,
        parsed: &'a [Option<Arc<ParsedSubject>>],
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
pub(crate) struct SubjectValues<'a> {
    id: Cow<'a, [u8]>,
    accession: Cow<'a, [u8]>,
    accession_version: Cow<'a, [u8]>,
    gi: i64,
    all_ids: Vec<DeflineIds>,
    titles: Option<Vec<&'a [u8]>>,
    length: usize,
    taxonomy: Cow<'a, TabularTaxonomy>,
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
        }
    }
}

impl<'a> TabularSubjects<'a> {
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
            TabularSubjects::Records {
                ids,
                length,
                parsed,
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
                    let [id, accession, accession_version] =
                        subject.tabular_ids(false).map_err(refused)?;
                    let all_ids = if needs.all_deflines {
                        subject.all_tabular_ids(false).map_err(refused)?
                    } else {
                        Vec::new()
                    };
                    return Ok(SubjectValues {
                        id: Cow::Owned(id),
                        accession: Cow::Owned(accession),
                        accession_version: Cow::Owned(accession_version),
                        gi: subject.tabular_gi(),
                        all_ids,
                        titles: None,
                        length,
                        taxonomy: Cow::Owned(TabularTaxonomy::default()),
                    });
                }
                let id: &[u8] = ids.get(s_idx).map_or(&b"unknown"[..], |id| id.as_ref());
                Ok(SubjectValues {
                    id: Cow::Borrowed(id),
                    accession: Cow::Borrowed(id),
                    accession_version: Cow::Borrowed(id),
                    gi: 0,
                    all_ids: Vec::new(),
                    titles: None,
                    length,
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
                                program,
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

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:844-851
// ```c++
//             	CRef<CBlast_def_line_set> bdlRef =
//             			CSeqDB::ExtractBlastDefline(subject_bh);
//             	if(setSubjectIds) {
//             		x_SetSubjectIds(subject_bh, bdlRef);
//             	}
//             	if(setSubjectTaxInfoAll) {
//             		x_SetTaxInfoAll(subject_bh, bdlRef);
//             	}
// ```
/// The taxonomy of the tabular fields (`staxid`, `staxids`;
/// `report::tabular::TabularTaxonomy::of_deflines`) of the reported sequences `ordinals` of `db`,
/// by rank: from the definition lines of a sequence with parsed Seq-ids (`parsed[rank]`), else
/// from its stored ones (`SeqDb::defline_set`, which keeps whether a line has a `taxid`), for
/// the ranks `ranks` (the subjects of the report's rows); the others are left empty. A header
/// that does not decode is refused before the report (`program` names the program).
pub(crate) fn tabular_taxonomy(
    db: &SeqDb,
    ordinals: &[u32],
    parsed: &[Option<Arc<ParsedSubject>>],
    ranks: impl IntoIterator<Item = usize>,
    program: &str,
) -> Result<Vec<TabularTaxonomy>> {
    let mut taxonomy = vec![TabularTaxonomy::default(); ordinals.len()];
    let mut read = vec![false; ordinals.len()];
    for rank in ranks {
        if rank >= ordinals.len() || std::mem::replace(&mut read[rank], true) {
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
    TabularWriter {
        spec,
        program_version,
        matrix: None,
    }
    .prolog(writer)?;
    Ok(())
}
