//! `PrintDomQuery` of a nucleotide query (rpstblastn's archives): the hits of each reading
//! frame with the frame after the query's label (`Query_1[-2]`), the sites with the codons of
//! their residues on the query's plus strand, and the motifs as ranges of bases.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1266-1280
//! ```c++
//! 	else	//na
//! 	{
//! 		(*m_ostream) << QUERY_TYPE_NUCLEOTIDE << DELIMIT << q.m_uiSeqLen << DELIMIT << q.m_strTitle << endl;
//!
//! 		vector<TDomSeqAlignIndex::__TCdAlignRecord> dimDomAligns[READINGFRAME::TOTAL_RFS], dimFeatVecs[READINGFRAME::TOTAL_RFS];
//!
//! 		bool bHasFeats = false, bHasRegFeats = false;
//!
//! 		for (READINGFRAME::TFRAMEINDEX rfidx = 0; rfidx < READINGFRAME::TOTAL_RFS; ++rfidx)
//! 		{
//! 			q.m_dimSplitAligns[rfidx].CreateRecordSets(q.m_vecAlignments, m_dom_cluster_src, dimDomAligns[rfidx], dimFeatVecs[rfidx], m_dmode);
//! 			bHasFeats = bHasFeats || !dimFeatVecs[rfidx].empty();
//! 			bHasRegFeats = (bHasRegFeats || !(dimFeatVecs[rfidx].empty() || dimFeatVecs[rfidx][0].pCdInfo->m_bIsStructDom));
//! 		}
//! ```

use std::collections::BTreeMap;

use super::args::TargetData;
use super::data::{self, CddData};
use super::domains::{
    domain_records, feature_records, map_seg_set, pr2na_ends, translated_pos_map, DataMode,
    DomAlign, FeatureRecord, Record,
};
use super::frames::{minus_pr_to_plus_na, plus_pr_to_plus_na, translate_all, RF_IDS, RF_SIZE};
use super::query::DomQuery;
use super::segset::SegSet;
use super::sites::MotifState;
use super::{
    mostly_mapped, print_cluster_line, print_dom_line, printed_position, site_segset, write_source,
    MotifRecord, PrintError, Report, OVERLAPLEADING,
};

/// The DOMAINS and SUPERFAMILIES parts of a nucleotide query (`-t doms` and `both`), then its
/// SITES and MOTIFS parts (`-t feats` and `both`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1281-1362
// ```c++
// 		if (TTargetData::e_feats != m_tdata)
// 		{
// 			map<int, TDomSeqAlignIndex::__TCdAlignRecord > dimDomFams[READINGFRAME::TOTAL_RFS];
//
// 			(*m_ostream) << DOMSTART << endl;
// 			for (READINGFRAME::TFRAMEINDEX rfidx = 0; rfidx < READINGFRAME::TOTAL_RFS; ++rfidx)
// 			{
// 				READINGFRAME::TFRAMEID rf = READINGFRAME::RF_IDS[rfidx];
// 				if (!dimDomAligns[rfidx].empty())
// 				{
// 					if (TDataModes::e_rep == m_dmode)
// 					{
// 						for (vector<TDomSeqAlignIndex::__TCdAlignRecord> :: const_iterator iterRec = dimDomAligns[rfidx].begin(), iterRecEnd = dimDomAligns[rfidx].end(); iterRecEnd != iterRec; ++iterRec)
// 						{
// 							if (NULL == iterRec->pCdInfo) continue;
// 							(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << '[' << rf << ']' << DELIMIT;
// 							if (1 == iterRec->pAlign->m_iRepClass)	//multi
// 								PrintDomLine(HITTYPE_MULTIDOM, *iterRec);
// 							else if (iterRec->pAlign->m_bSpecQualified)
// 							{
// 								PrintDomLine(HITTYPE_SPECIFIC, *iterRec);
// 								if (NULL != iterRec->pClst)
// 								{
// 									int key = iterRec->pAlign->m_iRegionIdx * OVERLAPLEADING + iterRec->pCdInfo->m_uiClusterPssmId;
// 									//if (dimDomFams[rfidx].end() == dimDomFams[rfidx].find(key))
// 									dimDomFams[rfidx].emplace(key, *iterRec);
// 								}
// 							}
// 							else if (NULL != iterRec->pClst)
// 								PrintClstLine(*iterRec);
// 							else	//not qualified for specific, but no cluster -- non-specific
// 								PrintDomLine(HITTYPE_NONSPECIFIC, *iterRec);
// 							(*m_ostream) << endl;
// 						}	//dom loops end
// 					}	//rep mode
// 					else	//nonrep mode
// 					{
// 						for (...)
// 						{
// 							if (NULL == iterRec->pCdInfo) continue;
// 							(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << '[' << rf << ']' << DELIMIT;
// 							if (1 == iterRec->pAlign->m_iRepClass)	//multi
// 								PrintDomLine(HITTYPE_MULTIDOM, *iterRec);
// 							else if (iterRec->pAlign->m_bRep && iterRec->pAlign->m_bSpecQualified)
// 								PrintDomLine(HITTYPE_SPECIFIC, *iterRec);
// 							else
// 								PrintDomLine(HITTYPE_NONSPECIFIC, *iterRec);
//
// 							if (NULL != iterRec->pClst)
// 							{
// 								int key = iterRec->pAlign->m_iRegionIdx * OVERLAPLEADING + iterRec->pCdInfo->m_uiClusterPssmId;
// 								//if (dimDomFams[rfidx].end() == dimDomFams[rfidx].find(key))
// 								dimDomFams[rfidx].emplace(key, *iterRec);
// 							}
//
// 							(*m_ostream) << endl;
// 						}
// 					}
// 				}
// 			}
// 			(*m_ostream) << DOMEND << endl;
//
// 			if (m_show_clusters)	//user want to see all superfams
// 			{
// 				(*m_ostream) << FAMSTART << endl;
// 				for (READINGFRAME::TFRAMEINDEX rfidx = 0; rfidx < READINGFRAME::TOTAL_RFS; ++rfidx)
// 				{
// 					//int rf = READINGFRAME::RF_IDS[rfidx];
//
// 					for (const auto & v : dimDomFams[rfidx])
// 					{
// 						(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << DELIMIT;
// 						PrintClstLine(v.second);
// 						(*m_ostream) << endl;
// 					}
//
// 				}
// 				(*m_ostream) << FAMEND << endl;
// 			}
// 		}
// ```
// Unlike a protein query's, the DOMAINS part (and SUPERFAMILIES with `-f`) is printed when no
// frame has a hit. The key of a superfamily is the `unsigned int` product of the region and
// `OVERLAPLEADING` (the cluster's PSSM-ID is 0, see `print_dom_query`) converted to `int`: the
// map orders the keys as signed numbers (region 2 comes before region 0). `Multidom` is never
// made (`domains::calculate`).
#[allow(clippy::too_many_arguments)]
pub(super) fn print_frames(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    mode: DataMode,
    show_families: bool,
    target: TargetData,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    if target != TargetData::Feats {
        report.line(&[b"DOMAINS"]);
        let mut families: [BTreeMap<i32, Record>; 6] = Default::default();
        for (rfidx, frame_families) in families.iter_mut().enumerate() {
            let rf = RF_IDS[rfidx].to_string();
            for record in &domain_records(&query.frames[rfidx], &query.aligns, data, mode) {
                let align = &query.aligns[record.align];
                report.write(session.as_bytes());
                report.write(b"\t");
                report.write(&query.accession);
                report.write(b"[");
                report.write(rf.as_bytes());
                report.write(b"]\t");
                let key = (align.region as u32).wrapping_mul(OVERLAPLEADING) as i32;
                match mode {
                    DataMode::Rep => {
                        if align.spec_qualified {
                            print_dom_line(report, b"Specific", record, query, data);
                            if record.cluster.is_some() {
                                frame_families.entry(key).or_insert(*record);
                            }
                        } else if record.cluster.is_some() {
                            print_cluster_line(report, record, query, data);
                        } else {
                            print_dom_line(report, b"Non-specific", record, query, data);
                        }
                    }
                    DataMode::Std | DataMode::Full => {
                        let kind: &[u8] = if align.rep && align.spec_qualified {
                            b"Specific"
                        } else {
                            b"Non-specific"
                        };
                        print_dom_line(report, kind, record, query, data);
                        if record.cluster.is_some() {
                            frame_families.entry(key).or_insert(*record);
                        }
                    }
                }
                report.write(b"\n");
            }
        }
        report.line(&[b"ENDDOMAINS"]);
        if show_families {
            report.line(&[b"SUPERFAMILIES"]);
            for record in families.iter().flat_map(|frame| frame.values()) {
                report.write(session.as_bytes());
                report.write(b"\t");
                report.write(&query.accession);
                report.write(b"\t");
                print_cluster_line(report, record, query, data);
                report.write(b"\n");
            }
            report.line(&[b"ENDSUPERFAMILIES"]);
        }
    }
    if target != TargetData::Doms {
        print_features(report, query, session, data, mode, motif_state)?;
    }
    Ok(())
}

/// The source column of a site line: the PSSM-ID or a model's accession (`write_source`), the
/// PSSM-ID as it is (the minus-strand frames), or `0` (the structural domains' sites).
#[derive(Clone, Copy)]
enum Source {
    Lookup(u32),
    Raw(u32),
    Zero,
}

/// The SITES and MOTIFS parts of a nucleotide query.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1364-1394
// ```c++
// 		if (TTargetData::e_doms != m_tdata)
// 		{
//
// 			if (bHasFeats)
// 			{
//
// 				const _TOfflDomQuery * p_ofl_q = dynamic_cast< const _TOfflDomQuery * > (&q);
// 				string __dimTranslated[READINGFRAME::TOTAL_RFS];
//
// 				const string *dimTranslated = __dimTranslated;
// 				if (nullptr == p_ofl_q)	//not _TOfflDomQuery
// 				{
// 					TranslateAll(q.m_strSeqData, q.m_iGenCode, __dimTranslated);
// 				}
// 				else
// 					dimTranslated = p_ofl_q->m_dimTranslated;
// 				...
// 				vector<TDomSeqAlignIndex::__MotifType_base> dimMotifTypes[READINGFRAME::TOTAL_RFS];
// 				...
// 				bool bHasMotifs = false;
//
// 				if (bHasRegFeats)
// 					(*m_ostream) << FEATSTART << endl;
//
// 				int iNegRF = READINGFRAME::TOTAL_RFS / 2;
// 				for (READINGFRAME::TFRAMEINDEX rfidx = 0; rfidx < iNegRF; ++rfidx)	//plus strand reading frames
// 				{
//
// 					READINGFRAME::TFRAMEID rf = READINGFRAME::Idx2Id(rfidx);
//
// 					vector<TDomSeqAlignIndex::__TCdAlignRecord> :: const_iterator iterFeatAlign = dimFeatVecs[rfidx].begin(), iterFeatAlignEnd = dimFeatVecs[rfidx].end();
// 					...
// 					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
// 					{
// 						if (iterFeatAlign->pCdInfo->m_bIsStructDom)
// 						{
// 							...
// 							break;
// 						}
//
// 						bool bIsSpecific = iterFeatAlign->pAlign->m_bRep && iterFeatAlign->pAlign->m_bSpecQualified;
// 						const char * pType = ANNOTTYPE_SPECIFIC;
// 						const list <TDomSite> * pFeats = &(iterFeatAlign->pCdInfo->m_lstSpecificFeatures);;
// 						unsigned int uiSrcPssmId = iterFeatAlign->pAlign->m_uiPssmId;
//
// 						if (!bIsSpecific && NULL != iterFeatAlign->pRootCdInfo)
// 						{
// 							pType = ANNOTTYPE_GENERIC;
// 							pFeats = &(iterFeatAlign->pCdInfo->m_lstGenericFeatures);
// 							uiSrcPssmId = iterFeatAlign->pCdInfo->m_uiHierarchyRoot;
// 						}
//
// 						for (list<TDomSite> :: const_iterator iterFeat = pFeats->begin(); iterFeat != pFeats->end(); ++iterFeat)
// 						{
// 							if (TDomSite::eType_StructMotif == iterFeat->m_iType)
// 							{
// 								dimMotifTypes[rfidx].emplace_back(iterFeat, iterFeatAlign, bIsSpecific, uiSrcPssmId);
// 								bHasMotifs = true;
// 								continue;
// 							}
// 							...
// 					// -- non-motif from sd domains
// 					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
// 					{
// 						...
// 						for (list<TDomSite> :: const_iterator iterFeat = iterFeatAlign->pCdInfo->m_lstSpecificFeatures.begin(); iterFeat != iterFeatAlign->pCdInfo->m_lstSpecificFeatures.end(); ++iterFeat)
// 						{
// 							if (TDomSite::eType_StructMotif == iterFeat->m_iType)
// 							{
// 								dimMotifTypes[rfidx].emplace_back(iterFeat, iterFeatAlign, true, 0);
// 								bHasMotifs = true;
// 								continue;
// 							}
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1573-1594,1700-1706
// ```c++
// 				for (READINGFRAME::TFRAMEINDEX rfidx = iNegRF; rfidx < READINGFRAME::TOTAL_RFS; ++rfidx)	//minus strand reading frames
// 				{
// 					READINGFRAME::TFRAMEID rf = READINGFRAME::Idx2Id(rfidx);
//
// 					vector<TDomSeqAlignIndex::__TCdAlignRecord> :: const_iterator iterFeatAlign = dimFeatVecs[rfidx].begin(), iterFeatAlignEnd = dimFeatVecs[rfidx].end();
// 					//SeqLen_t alnLen = iterFeatAlign->pAlign->m_ReadingFrame >> 2, alnLenAdjust = dimTranslated[rfidx].size() - alnLen;
//
// 					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
// 					{
// 						if (iterFeatAlign->pCdInfo->m_bIsStructDom)
// 						{
// 							...
// 							break;
// 						}
// 						... (as the plus strand's)
// 				}	//negative strand reading frames
//
// 				if (bHasRegFeats)
// 					(*m_ostream) << FEATEND << endl;
// ```
// The SITES header and footer depend on a regular (not structural) first record in some frame;
// the structural domains' sites of a frame are printed between them, or without them when no
// frame starts with a regular record. The translations are made only when some frame has a
// record (`bHasFeats`); `frames::translate_all` refuses where NCBI's `TranslateAll` erases from
// an empty string.
fn print_features(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    mode: DataMode,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let records: [Vec<FeatureRecord>; 6] = std::array::from_fn(|rfidx| {
        feature_records(&query.frames[rfidx], &query.aligns, data, mode)
    });
    if records.iter().all(Vec::is_empty) {
        return Ok(());
    }
    let has_regular = records.iter().any(|frame| {
        frame
            .first()
            .is_some_and(|record| !data.domains[record.domain].is_struct_dom)
    });
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1370-1379
    // ```c++
    // 				const _TOfflDomQuery * p_ofl_q = dynamic_cast< const _TOfflDomQuery * > (&q);
    // 				string __dimTranslated[READINGFRAME::TOTAL_RFS];
    //
    // 				const string *dimTranslated = __dimTranslated;
    // 				if (nullptr == p_ofl_q)	//not _TOfflDomQuery
    // 				{
    // 					TranslateAll(q.m_strSeqData, q.m_iGenCode, __dimTranslated);
    // 				}
    // 				else
    // 					dimTranslated = p_ofl_q->m_dimTranslated;
    // ```
    // (A query of the XML input is a `_TOfflDomQuery`: its frames' residues are those that its
    // alignments wrote, `offline::OfflineQuery`.)
    let translated;
    let translations: &[Vec<u8>; 6] = match &query.translated {
        Some(frames) => frames,
        None => {
            translated = translations(&query.residues)?;
            &translated
        }
    };
    if has_regular {
        report.line(&[b"SITES"]);
    }
    let mut motifs: [Vec<MotifRecord>; 6] = Default::default();
    for rfidx in 0..6 {
        let frame = Frame {
            rfidx: rfidx as u32,
            label: RF_IDS[rfidx].to_string(),
            translation: &translations[rfidx],
        };
        let mut at = 0usize;
        while at < records[rfidx].len() {
            let record = &records[rfidx][at];
            if data.domains[record.domain].is_struct_dom {
                break;
            }
            let align = &query.aligns[record.align];
            let is_specific = align.rep && align.spec_qualified;
            let (kind, specific, source): (&[u8], bool, u32) =
                if !is_specific && record.root.is_some() {
                    (
                        b"Generic",
                        false,
                        data.domains[record.domain].hierarchy_root,
                    )
                } else {
                    (b"Specific", true, align.pssm_id)
                };
            for site in data.sites(record.domain, specific) {
                if site.site_type == data::TYPE_STRUCT_MOTIF {
                    motifs[rfidx].push(MotifRecord {
                        site,
                        record: at,
                        source,
                    });
                    continue;
                }
                let source = if rfidx < 3 {
                    Source::Lookup(source)
                } else {
                    Source::Raw(source)
                };
                print_site(
                    report,
                    query,
                    session,
                    data,
                    &frame,
                    align,
                    site,
                    kind,
                    source,
                    motif_state,
                )?;
            }
            at += 1;
        }
        while at < records[rfidx].len() {
            let record = &records[rfidx][at];
            let align = &query.aligns[record.align];
            for site in data.sites(record.domain, true) {
                if site.site_type == data::TYPE_STRUCT_MOTIF {
                    motifs[rfidx].push(MotifRecord {
                        site,
                        record: at,
                        source: 0,
                    });
                    continue;
                }
                print_site(
                    report,
                    query,
                    session,
                    data,
                    &frame,
                    align,
                    site,
                    b"Specific",
                    Source::Zero,
                    motif_state,
                )?;
            }
            at += 1;
        }
    }
    if has_regular {
        report.line(&[b"ENDSITES"]);
    }
    if motifs.iter().any(|frame| !frame.is_empty()) {
        print_motifs(report, query, session, data, &records, &motifs)?;
    }
    Ok(())
}

/// `TranslateAll` of a query whose hits have sites or motifs, or the refusal where NCBI's
/// `TranslateAll` erases from an empty string (`frames::translate_all`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1376
// ```c++
// 					TranslateAll(q.m_strSeqData, q.m_iGenCode, __dimTranslated);
// ```
fn translations(residues: &[u8]) -> Result<[Vec<u8>; 6], PrintError> {
    translate_all(residues).ok_or_else(|| {
        PrintError::Unsupported(
            format!(
                "a nucleotide query of {} bases whose reverse complement has {} IUPACna letters, fewer than 3 in either (NCBI's TranslateAll erases from an empty string, undefined), is not supported by LOSAT's rpsbproc",
                residues.len(),
                super::frames::complement(residues).len()
            )
            .into_bytes(),
        )
    })
}

/// A reading frame of the site output: its index, its label and its translation.
struct Frame<'a> {
    rfidx: u32,
    label: String,
    translation: &'a [u8],
}

/// One site line of a nucleotide query: the residues of the frame's translation with the
/// plus-strand codon of each (`L33-35`; descending on the minus strand, `L35-33`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1435-1459
// ```c++
// 							CSegSet featsegs(*iterFeat);
// 							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
//
// 							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
// 							{
// 								vector<TSeg_base::TResiduePos> vecMappedPos;
// 								iterFeatAlign->pAlign->GetTranslatedPosMap(featsegs, q.m_uiSeqLen, vecMappedPos);
//
// 								if (iterFeat->MotifCheck(vecMappedPos, dimTranslated[rfidx]) > 0) continue;	//failed motif check
//
// 								(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << '[' << rf << ']' << DELIMIT << pType << DELIMIT << iterFeat->m_strTitle << DELIMIT;
//
// 								size_t mapped = 0;
// 								char dimDelimit[2] = {0, 0};
//
//
// 								for (vector<TSeg_base::TResiduePos> :: const_iterator iterRP = vecMappedPos.begin(), iterEnd = vecMappedPos.end(); iterEnd != iterRP; ++iterRP)
// 								{
// 									SeqPos_t xres = READINGFRAME::PlusPr2PlusNA(iterRP->curr, rfidx, q.m_uiSeqLen);
// 									(*m_ostream) << dimDelimit << (char)toupper(dimTranslated[rfidx][iterRP->curr]) << (xres + COORDSBASE) << '-' << (xres + codonEnd + COORDSBASE);
// 									dimDelimit[0] = COORDELIMIT;
// 								}
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1484-1496
// ```c++
// 								(*m_ostream) << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT;
// 								if (uiSrcPssmId < INVALIDPSSMID)
// 									(*m_ostream) << uiSrcPssmId;
// 								else
// 								{
// 									const TDomain *pSrcDom = m_dom_cluster_src.FindCdInfo(uiSrcPssmId);
// 									if (nullptr != pSrcDom)
// 										(*m_ostream) << pSrcDom->m_strAccession;
// 									else
// 										(*m_ostream) << uiSrcPssmId;
// 								}
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1628-1652
// ```c++
// 								for (vector<TSeg_base::TResiduePos> :: const_iterator iterRP = vecMappedPos.begin(), iterEnd = vecMappedPos.end(); iterEnd != iterRP; ++iterRP)
// 								{
// 									SeqPos_t xres = READINGFRAME::MinusPr2PlusNA(iterRP->curr, rfidx, q.m_uiSeqLen);
//
// 									(*m_ostream) << dimDelimit << (char)toupper(dimTranslated[rfidx][iterRP->curr]) << (xres + COORDSBASE) << '-' << (xres - codonEnd + COORDSBASE);
// 									dimDelimit[0] = COORDELIMIT;
// 								}
// 								...
// 								(*m_ostream) << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT << uiSrcPssmId << endl;
// ```
// `mapped` is never counted in this branch (0 in every line). `codonEnd` is a `SeqLen_t`
// (`unsigned int`), so the second number of a codon is an `unsigned int` (`xres + 2 + 1`,
// `xres - 2 + 1`, modulo 2^32) while the first is `xres + 1` in `int` (`printed_position`). The
// structural domains' sites (`rpsbproc_main.cpp:1499-1569`, `:1658-1699`)
// print `Specific` and the source `0`; the minus-strand frames print the source PSSM-ID as it
// is, the plus-strand frames as `write_source`. A residue position outside the frame's
// translation other than its end (where `operator[]` gives the terminating 0) reads outside
// the string (undefined): refused.
#[allow(clippy::too_many_arguments)]
fn print_site(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    frame: &Frame,
    align: &DomAlign,
    site: &data::Site,
    kind: &[u8],
    source: Source,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let mut featsegs = site_segset(data, site);
    map_seg_set(&mut featsegs, align);
    if !mostly_mapped(&featsegs, site) {
        return Ok(());
    }
    let mapped_positions = translated_pos_map(&featsegs, align, query.length)
        .map_err(|what| PrintError::Exception(what.into_bytes()))?;
    if motif_state
        .motif_check(data, site, &mapped_positions, frame.translation)
        .map_err(|error| super::motif_error(error, data, site))?
        > 0
    {
        return Ok(());
    }
    let mut codons = Vec::with_capacity(mapped_positions.len() * 12);
    let mut delimit: &[u8] = b"";
    for position in &mapped_positions {
        let letter = translated_residue(frame, position.curr)?;
        let (first, last) = codon(frame, position.curr, query.length);
        codons.extend_from_slice(delimit);
        codons.push(letter);
        codons.extend_from_slice(format!("{first}-{last}").as_bytes());
        delimit = b",";
    }
    report.write(session.as_bytes());
    report.write(b"\t");
    report.write(&query.accession);
    report.write(b"[");
    report.write(frame.label.as_bytes());
    report.write(b"]\t");
    report.write(kind);
    report.write(b"\t");
    report.write(data.text(site.title));
    report.write(b"\t");
    report.write(&codons);
    let complete: i64 = data
        .site_segs(site)
        .iter()
        .map(|seg| i64::from(seg[1]) - i64::from(seg[0]) + 1)
        .sum();
    report.write(format!("\t{complete}\t0\t").as_bytes());
    match source {
        Source::Lookup(source) => write_source(report, data, source),
        Source::Raw(source) => report.write(source.to_string().as_bytes()),
        Source::Zero => report.write(b"0"),
    }
    report.write(b"\n");
    Ok(())
}

/// `(char)toupper(dimTranslated[rfidx][curr])`: the terminating 0 at the translation's end;
/// elsewhere outside it NCBI reads outside the string (undefined), which LOSAT refuses.
fn translated_residue(frame: &Frame, curr: i64) -> Result<u8, PrintError> {
    let len = frame.translation.len();
    match usize::try_from(curr) {
        Ok(at) if at < len => Ok(frame.translation[at].to_ascii_uppercase()),
        Ok(at) if at == len => Ok(0),
        _ => Err(PrintError::Unsupported(
            format!(
                "a site at residue {} of the translation of reading frame {} ({len} residues; NCBI reads outside it) is not supported by LOSAT's rpsbproc",
                printed_position(curr),
                frame.label
            )
            .into_bytes(),
        )),
    }
}

/// The two numbers of the codon of residue `curr` of the frame's translation (quoted at
/// `print_site`): `xres + COORDSBASE` in `int`, `xres + codonEnd + COORDSBASE` (minus strand
/// `xres - codonEnd + COORDSBASE`) in `unsigned int`. The `int` arithmetic that overflows (`pr -
/// 1` of `MinusPr2PlusNA` at `INT_MIN`, `xres + 1` at `INT_MAX`; undefined) is computed as its
/// value (`frames::minus_pr_to_plus_na`, `printed_position`; PD-LOSAT-NCBI-DEFECTS exception 11).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:949,1456-1457,1631-1633
// ```c++
// 	static constexpr const SeqLen_t codonEnd = READINGFRAME::RF_SIZE - 1;
// 									SeqPos_t xres = READINGFRAME::PlusPr2PlusNA(iterRP->curr, rfidx, q.m_uiSeqLen);
// 									(*m_ostream) << dimDelimit << (char)toupper(dimTranslated[rfidx][iterRP->curr]) << (xres + COORDSBASE) << '-' << (xres + codonEnd + COORDSBASE);
// 									SeqPos_t xres = READINGFRAME::MinusPr2PlusNA(iterRP->curr, rfidx, q.m_uiSeqLen);
//
// 									(*m_ostream) << dimDelimit << (char)toupper(dimTranslated[rfidx][iterRP->curr]) << (xres + COORDSBASE) << '-' << (xres - codonEnd + COORDSBASE);
// ```
fn codon(frame: &Frame, curr: i64, query_len: u32) -> (i64, u32) {
    let (xres, last) = if frame.rfidx < RF_SIZE {
        let xres = plus_pr_to_plus_na(curr, frame.rfidx);
        (xres, (xres as u32).wrapping_add(2).wrapping_add(1))
    } else {
        let xres = minus_pr_to_plus_na(curr, frame.rfidx, query_len);
        (xres, (xres as u32).wrapping_sub(2).wrapping_add(1))
    };
    (printed_position(xres), last)
}

/// The MOTIFS part of a nucleotide query: each motif as the range of bases of its mapped
/// residues.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1709-1795
// ```c++
// 				if (bHasMotifs)
// 				{
// 					(*m_ostream) << MOTIFSTART << endl;
//
// 					// -- plus strand
// 					for (READINGFRAME::TFRAMEINDEX rfidx = 0; rfidx < iNegRF; ++rfidx)
// 					{
// 						READINGFRAME::TFRAMEID rf = READINGFRAME::Idx2Id(rfidx);
// 						for (vector<TDomSeqAlignIndex::__MotifType_base> :: const_iterator iterM = dimMotifTypes[rfidx].begin(), iterMEnd = dimMotifTypes[rfidx].end(); iterMEnd != iterM; ++iterM)
// 						{
// 							CSegSet featsegs(*iterM->iterMotifFeat);
// 							iterM->iterAlignRec->pAlign->MapSegSet(featsegs, false);
//
// 							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterM->iterMotifFeat->GetCompleteSize() > 0.8))
// 							{
// 								iterM->iterAlignRec->pAlign->Pr2NaConvert(featsegs);
// 								(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << '[' << rf << ']' << DELIMIT << iterM->iterMotifFeat->m_strTitle << DELIMIT << featsegs.GetLeft() + COORDSBASE << DELIMIT << featsegs.GetRight() + codonEnd + COORDSBASE << DELIMIT << iterM->uiSrcPssmId << endl;
// 								...
// 							}
// 						}
// 						...
// 					}
//
// 					// -- minus strand
// 					for (READINGFRAME::TFRAMEINDEX rfidx = iNegRF; rfidx < READINGFRAME::TOTAL_RFS; ++rfidx)
// 					{
// 						int rf = 2 - rfidx;
// 						vector<TDomSeqAlignIndex::__TCdAlignRecord> :: const_iterator iterFeatAlign = dimFeatVecs[rfidx].begin(), iterFeatAlignEnd = dimFeatVecs[rfidx].end();
// 						SeqLen_t alnLen = iterFeatAlign->pAlign->m_ReadingFrame >> 2, alnLenAdjust = dimTranslated[rfidx].size() - alnLen;
//
// 						for (vector<TDomSeqAlignIndex::__MotifType_base> :: const_iterator iterM = dimMotifTypes[rfidx].begin(), iterMEnd = dimMotifTypes[rfidx].end(); iterMEnd != iterM; ++iterM)
// 						{
// 							CSegSet featsegs(*iterM->iterMotifFeat);
// 							iterM->iterAlignRec->pAlign->MapSegSet(featsegs, false);
//
// 							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterM->iterMotifFeat->GetCompleteSize() > 0.8))
// 							{
// 								iterM->iterAlignRec->pAlign->Pr2NaConvert(featsegs);
// 								(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << '[' << rf << ']' << DELIMIT << iterM->iterMotifFeat->m_strTitle << DELIMIT << featsegs.GetLeft() + COORDSBASE << DELIMIT << featsegs.GetRight() - codonEnd + COORDSBASE << DELIMIT << iterM->uiSrcPssmId << endl;
//
// 							}
// 						}
//
// 					}
// 					(*m_ostream) << MOTIFEND << endl;
// 				}
// ```
// No residue letter and no clipping of the structural domains' motifs (unlike a protein
// query's); the source is the PSSM-ID as it is. The minus-strand loop reads `alnLen` from the
// first record of a frame that may have none (survey 8.7): unused, so dead (no motif of that
// frame exists then), and not read here. `GetRight() +/- codonEnd + COORDSBASE` is an
// `unsigned int`; `GetLeft() + COORDSBASE` an `int` (`printed_position`).
fn print_motifs(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    records: &[Vec<FeatureRecord>; 6],
    motifs: &[Vec<MotifRecord>; 6],
) -> Result<(), PrintError> {
    report.line(&[b"MOTIFS"]);
    for (rfidx, frame_motifs) in motifs.iter().enumerate() {
        let rf = RF_IDS[rfidx].to_string();
        for motif in frame_motifs {
            let align = &query.aligns[records[rfidx][motif.record].align];
            let mut featsegs = site_segset(data, motif.site);
            map_seg_set(&mut featsegs, align);
            if !mostly_mapped(&featsegs, motif.site) {
                continue;
            }
            let (first, last) = motif_ends(&featsegs, align, rfidx);
            report.write(session.as_bytes());
            report.write(b"\t");
            report.write(&query.accession);
            report.write(b"[");
            report.write(rf.as_bytes());
            report.write(b"]\t");
            report.write(data.text(motif.site.title));
            report.write(format!("\t{first}\t{last}\t{}\n", motif.source).as_bytes());
        }
    }
    report.line(&[b"ENDMOTIFS"]);
    Ok(())
}

/// The two numbers of a motif line of frame `rfidx` (quoted at `print_motifs`):
/// `GetLeft() + COORDSBASE` in `int`, `GetRight() + codonEnd + COORDSBASE` (minus strand
/// `- codonEnd`) in `unsigned int`, after `Pr2NaConvert` of the set (not empty: `mostly_mapped`).
/// The `int` arithmetic that overflows (`domains::pr2na_ends`,
/// and `+ 1` at `INT_MAX`; undefined) is computed as its value (PD-LOSAT-NCBI-DEFECTS exception 11).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1746-1747,1788-1789
// ```c++
// 								iterM->iterAlignRec->pAlign->Pr2NaConvert(featsegs);
// 								(*m_ostream) << ... << featsegs.GetLeft() + COORDSBASE << DELIMIT << featsegs.GetRight() + codonEnd + COORDSBASE << DELIMIT << iterM->uiSrcPssmId << endl;
// 								iterM->iterAlignRec->pAlign->Pr2NaConvert(featsegs);
// 								(*m_ostream) << ... << featsegs.GetLeft() + COORDSBASE << DELIMIT << featsegs.GetRight() - codonEnd + COORDSBASE << DELIMIT << iterM->uiSrcPssmId << endl;
// ```
fn motif_ends(featsegs: &SegSet, align: &DomAlign, rfidx: usize) -> (i64, u32) {
    let (left, right) = pr2na_ends(featsegs, align).expect("a mapped motif is not empty");
    let last = if rfidx < 3 {
        (right as u32).wrapping_add(2).wrapping_add(1)
    } else {
        (right as u32).wrapping_sub(2).wrapping_add(1)
    };
    (printed_position(left), last)
}

#[cfg(test)]
#[path = "nucleotide_tests.rs"]
mod tests;
