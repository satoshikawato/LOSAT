//! The definition lines of a sequence read with `-parse_seqids`, as CWriteDB builds them
//! (`CWriteDB_Impl::ExtractBioseqDeflines`, `c++/src/objtools/blast/seqdb_writer/
//! writedb_impl.cpp`): the second parse of the raw definition line, split at ^A
//! (`x_GetFastaReaderDeflines`), the choice between it and the Bioseq's ids
//! (`s_UseFastaReaderDeflines`, `x_BuildDeflinesFromBioseq`), and the edits of
//! `CBuildDatabase` around it (the ids of the other molecule, `FixTaxId`).

use super::bioseq::{Bioseq, ParsedIds};
use super::fasta_source::{fasta_reader_fasta, post_seq_id_warning, Warn};
use crate::seqid::accguide::{identify_seq_id, parse_flags, F_NUC, F_PROT};
use crate::seqid::label::content;
use crate::seqid::model::ParsedSeqId;
use crate::seqid::parse::{parse_ids, seq_id_from_text, set_typed, SeqIdError};
use crate::seqid::rank::{best_rank_score, find_best};
use crate::seqid::{BlastDefline, DeflineSet};

/// `GetSeqIdString(true)`.
fn seq_id_string_with_version(id: &ParsedSeqId) -> Vec<u8> {
    let mut out = Vec::new();
    content(id, true, &mut out);
    out
}

/// `CSeq_id(CSeq_id::e_Local, label)`: the typed `Set` of a local id (`SetStrOrId` of the
/// trimmed label), which does not fail.
fn local_from_label(label: &[u8]) -> ParsedSeqId {
    set_typed(crate::seqid::accguide::SEQID_LOCAL, label, b"", 0, b"")
        .ok()
        .flatten()
        .expect("a local id is always made")
}

/// Whether `IdentifyAccession()` of `id` names the other molecule than `is_protein`:
/// `(info & (fAcc_prot | fAcc_nuc)) && is_protein == !!(info & fAcc_nuc)`.
fn of_other_molecule(id: &ParsedSeqId, is_protein: bool, warn: Warn<'_>) -> bool {
    let info = identify_seq_id(
        id,
        parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
        &mut |w| post_seq_id_warning(w, &mut *warn),
    );
    info & (F_PROT | F_NUC) != 0 && is_protein == (info & F_NUC != 0)
}

/// The ids edit of `CBuildDatabase::AddSequences` (without `[BLAST] LONG_SEQID`): every
/// id that is neither local nor a GI and whose accession is of the other molecule
/// becomes the local id of its label.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:842-856
/// ```c++
///             if (!m_LongIDs) {
///
///                 // If accession's molecule type is different than expected,
///                 // change sequence id to local. CFastaReader cannot distingush
///                 // between bare pir protein ids genbank nucleotide ids.
///                 CBioseq* bss = const_cast<CBioseq*>(bs.GetNonNullPointer());
///                 for (auto& it: bss->SetId()) {
///                     CSeq_id::EAccessionInfo info = it->IdentifyAccession();
///                     if (!it->IsLocal() && !it->IsGi() &&
///                         (info & (CSeq_id::fAcc_prot | CSeq_id::fAcc_nuc)) &&
///                         m_IsProtein == !!(info & CSeq_id::fAcc_nuc)) {
///
///                         string label = it->GetSeqIdString(true);
///                         it.Reset(new CSeq_id(CSeq_id::e_Local, label));
///                     }
///                 }
///             }
/// ```
/// `IdentifyAccession()` runs for every id, local and GI ones included (its diagnostic,
/// if any, is posted for them too).
pub fn localize_other_molecule(ids: &mut [ParsedSeqId], is_protein: bool, warn: Warn<'_>) {
    for id in ids {
        let other = of_other_molecule(id, is_protein, warn);
        if !matches!(id, ParsedSeqId::Local(_) | ParsedSeqId::Gi(_)) && other {
            *id = local_from_label(&seq_id_string_with_version(id));
        }
    }
}

/// The GI check of `CBuildDatabase::AddSequences`: when an id is a GI above 4294967295,
/// the sequence is skipped and the GI of the last GI id is logged.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:812-833
/// ```c++
/// 	    if ( m_SkipLargeGis && !ids.empty() && ids.front().NotEmpty()){
/// 		bool skip_this = false;
/// 		for(list< CRef<CSeq_id> >::const_iterator it = ids.begin(); it != ids.end(); it++  ){
/// 		    if( it->NotEmpty() ){
/// 			CSeq_id::EAccessionInfo info = (*it)->IdentifyAccession();
/// 			if( info == CSeq_id::EAccessionInfo::eAcc_gi ){
/// 			    check_gi = (*it)->GetGi();
/// 			    if( check_gi > max_gi32_val )
/// 			    {
/// 				skip_this = true;
/// 			    }
/// 			}
/// 		    }
/// 		}
/// 		if( skip_this ){
/// 		    m_LogFile << "Ignoring gi '" << check_gi << "' as it has value larger then " << 0xFFFFFFFF<< endl;
/// 		    bs = src.GetNext();
/// 		    continue;
/// 		}
/// 	    }
/// ```
/// `IdentifyAccession()` of a GI is `eAcc_gi` (its `E_Choice`), of no other id. The
/// diagnostics of `IdentifyAccession` for the other ids are posted to `warn`.
pub fn large_gi_to_skip(ids: &[ParsedSeqId], warn: Warn<'_>) -> Option<i64> {
    let mut skip_this = false;
    let mut check_gi = 0;
    for id in ids {
        identify_seq_id(
            id,
            parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
            &mut |w| post_seq_id_warning(w, &mut *warn),
        );
        if let ParsedSeqId::Gi(gi) = *id {
            check_gi = gi;
            if gi > 0xFFFF_FFFF {
                skip_this = true;
            }
        }
    }
    skip_this.then_some(check_gi)
}

/// The second parse of `x_GetFastaReaderDeflines(bioseq, ..., accept_gt = false,
/// parse_ids = true, long_seqids = false, scan = true)`: one definition line per part of
/// the raw definition line between ^A; `None` where NCBI's parse throws a
/// `CSeqIdException` (which `x_ExtractDeflines` catches) or finds no definition line.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1766-1767,1807-1880
/// ```c++
///     if (fasta.empty())
///         return;
///     ...
///     } else {
///
///         int skip = 1;
///         while(fasta.size()) {
///             size_t id_start = skip;
///             size_t pos_title = fasta.find(" ", skip);
///             size_t pos_next = fasta.find("\001", skip);
///             skip = 1;
///
///             if (pos_next == fasta.npos) {
///                 if (accept_gt) {
///                     pos_next = fasta.find(" >");
///                     skip = 2;
///                 }
///             } else {
///                 // If there is a ^A, turn off GT checking.
///                 accept_gt = false;
///             }
///
///             if (pos_next == fasta.npos) {
///                 pos_next = fasta.size();
///                 skip = 0;
///             }
///
///             if (pos_title == fasta.npos || pos_title >= pos_next) {
///                 // title field is missing
///                 pos_title = pos_next;
///             }
///
///             string ids(fasta, id_start, pos_title - id_start);
///             if (pos_title == pos_next) pos_title--;
///             string title(fasta, pos_title + 1, pos_next-pos_title - 1);
///             string remaining(fasta, pos_next, fasta.size() - pos_next);
///             fasta.swap(remaining);
///
///             // Parse ids.  They may or may not be bar-separated.
///             list< CRef<CSeq_id> > seqids;
///             if (ids.find('|') != NPOS){
///             	CSeq_id::ParseFastaIds(seqids, ids);
///             }
///             else {
///             	CRef<CSeq_id> id(new CSeq_id(ids, CSeq_id::fParse_RawText | CSeq_id::fParse_ValidLocal));
///             	if ((id->Which() == CSeq_id::e_Prf) ||
///             	    (id->Which() == CSeq_id::e_Pir)){
///             		string label = id->GetSeqIdString(true);
///                     id.Reset(new CSeq_id(CSeq_id::e_Local, label));
///                 }
///
///             	CSeq_id::EAccessionInfo info = id->IdentifyAccession();
///             	if ((info & (CSeq_id::fAcc_prot | CSeq_id::fAcc_nuc)) &&
///             	    (bioseq.IsAa() == !!(info & CSeq_id::fAcc_nuc))) {
///             	     string label = id->GetSeqIdString(true);
///             	     id.Reset(new CSeq_id(CSeq_id::e_Local, label));
///             	}
///
///             	seqids.push_back(id);
///             }
///
///             // Build the actual defline.
///
///             defline.Reset(new CBlast_def_line);
///             defline->SetSeqid().swap(seqids);
///             defline->SetTitle(title);
///             ...
///             bdls->Set().push_back(defline);
///         }
///     }
///     s_CheckEmptyLists(bdls, true);
///     deflines = bdls;
/// ```
/// `ParseFastaIds(seqids, ids)` is `ParseIDs` with `fParse_RawText | fParse_AnyLocal`
/// (no `fParse_PartialOK`, Seq_id.cpp:2603-2611): any exception goes up. No membership,
/// link or PIG is given (`ExtractBioseqDeflines` passes none). The `fasta` string is the
/// one `fasta_reader_title` describes; its escapes do not fail.
pub fn second_parse(parsed: &ParsedIds, is_aa: bool, warn: Warn<'_>) -> Option<DeflineSet> {
    let mut fasta = fasta_reader_fasta(&parsed.defline_field).ok()?;
    if fasta.is_empty() {
        return None;
    }
    let find = |text: &[u8], byte: u8, from: usize| {
        text.get(from..)
            .and_then(|rest| rest.iter().position(|&c| c == byte))
            .map(|offset| from + offset)
    };
    let mut deflines = Vec::new();
    let mut skip = 1usize;
    while !fasta.is_empty() {
        let id_start = skip;
        let mut pos_title = find(&fasta, b' ', skip);
        let mut pos_next = find(&fasta, 0x01, skip);
        skip = 1;
        if pos_next.is_none() {
            pos_next = Some(fasta.len());
            skip = 0;
        }
        let pos_next = pos_next.expect("set above");
        if pos_title.is_none_or(|pos| pos >= pos_next) {
            pos_title = Some(pos_next);
        }
        let mut pos_title = pos_title.expect("set above");
        let ids = fasta[id_start.min(pos_title)..pos_title].to_vec();
        if pos_title == pos_next {
            // `pos_title--` of 0 wraps; the title is then from 0 (`pos_title + 1`).
            pos_title = pos_title.wrapping_sub(1);
        }
        let title_start = pos_title.wrapping_add(1);
        let title = fasta[title_start.min(pos_next)..pos_next].to_vec();
        fasta.drain(..pos_next);
        let seqids = match parse_second_ids(&ids, is_aa, warn) {
            Ok(seqids) => seqids,
            Err(_) => return None,
        };
        deflines.push(BlastDefline {
            title: Some(title),
            ids: seqids,
            ..BlastDefline::default()
        });
    }
    Some(DeflineSet { deflines })
}

/// The ids of one part of the second parse (`x_GetFastaReaderDeflines`, quoted there).
fn parse_second_ids(
    ids: &[u8],
    is_aa: bool,
    warn: Warn<'_>,
) -> Result<Vec<ParsedSeqId>, SeqIdError> {
    let mut on_warning = |w| post_seq_id_warning(w, &mut *warn);
    if ids.contains(&b'|') {
        let parsed = parse_ids(
            ids,
            parse_flags::RAW_TEXT | parse_flags::ANY_LOCAL,
            &mut on_warning,
        )?;
        // Without fParse_PartialOK no id stays unset (an exception goes up first).
        return Ok(parsed.into_iter().flatten().collect());
    }
    let mut id = seq_id_from_text(
        ids,
        parse_flags::RAW_TEXT | parse_flags::VALID_LOCAL,
        &mut on_warning,
    )?;
    if matches!(id, ParsedSeqId::Prf(_) | ParsedSeqId::Pir(_)) {
        id = local_from_label(&seq_id_string_with_version(&id));
    }
    if of_other_molecule(&id, is_aa, warn) {
        id = local_from_label(&seq_id_string_with_version(&id));
    }
    Ok(vec![id])
}

/// `CBioseq::GetNonLocalId()` of a FASTA Bioseq: its best id by `s_BestNonLocalRank` if
/// that is not local, else the non-local row of its Seq-hist alignment (an id with a
/// range).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Bioseq.cpp:280-323
/// ```c++
/// static int s_BestNonLocalRank(const CRef<CSeq_id>& id)
/// {
///     if (id.Empty()) {
///         return kMax_Int;
///     } else if (id->IsLocal()) {
///         return kMax_Int - 1;
///     } else {
///         return id->BestRankScore();
///     }
/// }
///
/// const CSeq_id* CBioseq::GetNonLocalId() const
/// {
///     CRef<CSeq_id> id = FindBestChoice(GetId(), &s_BestNonLocalRank);
///     if (id.Empty()) {
///         return NULL; // No way to verify potential IDs found elsewhere
///     } else if ( !id->IsLocal() ) {
///         return &*id;
///     }
///
///     const CSeq_inst& inst = GetInst();
///     if ( !inst.CanGetHist()  ||  !inst.GetHist().CanGetAssembly() ) {
///         return NULL;
///     }
///
///     ITERATE (CSeq_hist::TAssembly, it, inst.GetHist().GetAssembly() ) {
///         try {
///             if ((*it)->CheckNumRows() != 2) {
///                 continue;
///             }
///         } catch (CSeqalignException&) { // fails basic validation; ignore
///             continue;
///         }
///
///         const CSeq_id& id1 = (*it)->GetSeq_id(0);
///         const CSeq_id& id2 = (*it)->GetSeq_id(1);
///         if (id1.IsLocal()  &&  id1.Match(*id)  &&  !id2.IsLocal()) {
///             return &id2;
///         } else if (id2.IsLocal()  &&  id2.Match(*id)  &&  !id1.IsLocal()) {
///             return &id1;
///         }
///     }
///
///     return NULL;
/// }
/// ```
/// The alignment of `CFastaReader::xCreateAlignment` has the generated `lcl|N` as its
/// first row and the replaced id as its second; the generated id is the Bioseq's only id
/// (the edits of `AddSequences` keep local ids), so the second row is returned when it is
/// not local.
fn non_local_id(parsed: &ParsedIds) -> Option<&ParsedSeqId> {
    let score = |id: &ParsedSeqId| {
        if matches!(id, ParsedSeqId::Local(_)) {
            i32::MAX - 1
        } else {
            best_rank_score(id)
        }
    };
    let best = find_best(&parsed.ids, score)?;
    if !matches!(best, ParsedSeqId::Local(_)) {
        return Some(best);
    }
    let old = parsed.hist_old_id.as_ref()?;
    let generated = parsed.ids.first()?;
    (generated == best && !matches!(old, ParsedSeqId::Local(_))).then_some(old)
}

/// `s_UseFastaReaderDeflines(bioseq, deflines, long_seqid = false)`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:736-762
/// ```c++
/// static bool s_UseFastaReaderDeflines(CConstRef<CBioseq> & bioseq, CConstRef<CBlast_def_line_set> & deflines, bool long_seqid)
/// {
/// 	if(deflines.Empty())
/// 		return false;
///
/// 	const CSeq_id * bioseq_id = bioseq->GetNonLocalId();
///
/// 	if(bioseq_id == NULL ||
///        // For bare pir and prf ids go with the one from defline.
///        // This is to parse bare ids as local ones. The bare pdb ids are pdb in
///        // bioseq (parsed by CFastaReader), but local in deflines (parsed by
///        // CSeq_id).
///        (!long_seqid && (bioseq_id->IsPrf() || bioseq_id->IsPir()))) {
///         return true;
///        }
///
/// 	// Bioseq has non-local id, make sure at least one id is non-local from CFastaReader
/// 	// defline
///     ITERATE(list< CRef<CBlast_def_line> >, iter, deflines->Get()) {
///         CRef<CSeq_id> id = FindBestChoice((**iter).GetSeqid(), &CSeq_id::BestRank);
///         if (id.NotEmpty()  &&  !id->IsLocal()) {
///                 return true;
///         }
///     }
///     return false;
///
/// }
/// ```
fn use_fasta_reader_deflines(parsed: &ParsedIds, deflines: &DeflineSet) -> bool {
    match non_local_id(parsed) {
        None | Some(ParsedSeqId::Prf(_)) | Some(ParsedSeqId::Pir(_)) => true,
        Some(_) => deflines.deflines.iter().any(|defline| {
            find_best(&defline.ids, best_rank_score)
                .is_some_and(|id| !matches!(id, ParsedSeqId::Local(_)))
        }),
    }
}

/// `CWriteDB::ExtractBioseqDeflines(bioseq, parse_ids = true, long_ids = false, scan =
/// true)` for a Bioseq of the FASTA path (after the edits of `AddSequences`), without the
/// serialization (the writer posts its diagnostics).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:802-836
/// ```c++
///         if (bin_hdr.empty()) {
/// 	    try {
///             x_GetFastaReaderDeflines(*bioseq,
///                                      deflines,
///                                      membbits,
///                                      linkouts,
///                                      pig,
///                                      false,
///                                      parse_ids,
///                                      long_ids,
///                                      scan_bioseq_4_cfastareader_usrobj);
/// 	      } catch ( const CSeqIdException&e ) {
/// 		//LOG_POST(Info << "x_GetFastaReaderDeflines " << e.GetMsg() );
/// 	    }
///         }
///
///         if(!s_UseFastaReaderDeflines(bioseq, deflines, long_ids)) {
///         	deflines.Reset();
///         }
///
///         if (bin_hdr.empty() && deflines.Empty()) {
///             x_BuildDeflinesFromBioseq(*bioseq,
///                                       deflines,
///                                       membbits,
///                                       linkouts,
///                                       pig);
///         }
///     }
///
///     if (bin_hdr.empty() &&
///         (deflines.Empty() || deflines->Get().empty())) {
///
///         NCBI_THROW(CWriteDBException,
///                    eArgErr,
///                    "Error: No deflines provided.");
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:606-614,650-656,694-695,716-719
/// ```c++
///     ITERATE(list< CRef< CSeqdesc > >, iter, bioseq.GetDescr().Get()) {
///         const CSeqdesc & desc = **iter;
///
///         if (desc.IsTitle()) {
///             //defline->SetTitle((**iter)->GetTitle());
///             titles = (**iter).GetTitle();
///         }
///     ...
///     list< CRef<CSeq_id> > ids = bioseq.GetId();
///     ...
///     while(! ids.empty()) {
///         defline.Reset(new CBlast_def_line);
///
///         defline->SetSeqid() = ids;
///         ids.clear();
///         ...
///         defline->SetTitle(titles);
///         ...
///         bdls->Set().push_back(defline);
///     }
///
///     s_CheckEmptyLists(bdls, true);
///     deflines = bdls;
/// ```
/// (`x_BuildDeflinesFromBioseq`: the one definition line of the Bioseq's ids and title
/// descriptor; a FASTA Bioseq has no taxon.) The Bioseq always has an id, so the "No
/// deflines" error does not occur.
pub fn extract_bioseq_deflines(bioseq: &Bioseq, warn: Warn<'_>) -> DeflineSet {
    let parsed = bioseq
        .parsed
        .as_ref()
        .expect("a Bioseq read with -parse_seqids");
    if let Some(deflines) = second_parse(parsed, bioseq.is_aa, warn) {
        if use_fasta_reader_deflines(parsed, &deflines) {
            return deflines;
        }
    }
    DeflineSet {
        deflines: vec![BlastDefline {
            title: Some(parsed.title_descr.clone().unwrap_or_default()),
            ids: parsed.ids.clone(),
            ..BlastDefline::default()
        }],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqid::model::{ObjectId, TextseqId};

    fn parsed(defline: &[u8], ids: Vec<ParsedSeqId>, title: Option<&[u8]>) -> ParsedIds {
        // The DefLine field is PrintableString(defline); for these lines (no backslash,
        // no control byte but ^A and tab) that is the line with \1 and \t escapes.
        let mut field = Vec::new();
        for &c in defline {
            match c {
                0x01 => field.extend_from_slice(b"\\1"),
                b'\t' => field.extend_from_slice(b"\\t"),
                _ => field.push(c),
            }
        }
        ParsedIds {
            ids,
            hist_old_id: None,
            defline_field: field,
            title_descr: title.map(<[u8]>::to_vec),
        }
    }

    fn refseq(acc: &str, version: i32) -> ParsedSeqId {
        ParsedSeqId::Other(TextseqId {
            accession: Some(acc.as_bytes().to_vec()),
            version: Some(version),
            ..TextseqId::default()
        })
    }

    fn titles_and_ids(set: &DeflineSet) -> Vec<(Vec<u8>, Vec<ParsedSeqId>)> {
        set.deflines
            .iter()
            .map(|d| (d.title.clone().unwrap(), d.ids.clone()))
            .collect()
    }

    // survey B fixture p02_multi (makeblastdb 2.17.0, blastdbcmd %i|%t).
    #[test]
    fn the_second_parse_splits_at_ctrl_a_and_keeps_the_spaces_of_titles() {
        let mut no_warn = |_: &[u8]| panic!("no warning");
        let p = parsed(
            b">ref|NP_100003.1| first title\x01ref|NP_100004.1| second title\x01gi|1005| third title",
            vec![refseq("NP_100003", 1)],
            Some(b"first title"),
        );
        let set = second_parse(&p, true, &mut no_warn).unwrap();
        assert_eq!(
            titles_and_ids(&set),
            [
                (b"first title".to_vec(), vec![refseq("NP_100003", 1)]),
                (b"second title".to_vec(), vec![refseq("NP_100004", 1)]),
                (b"third title".to_vec(), vec![ParsedSeqId::Gi(1005)]),
            ]
        );
        let p = parsed(
            b">ref|NP_100011.1|\ttab right after id",
            vec![refseq("NP_100011", 1)],
            None,
        );
        let set = second_parse(&p, true, &mut no_warn).unwrap();
        assert_eq!(titles_and_ids(&set)[0].0, b"  tab right after id");
        // A space right after '>': the id is empty, NCBI's parse throws.
        let p = parsed(
            b">  ref|NP_100008.1|  leading spaces",
            vec![refseq("NP_100008", 1)],
            Some(b"leading spaces"),
        );
        assert!(second_parse(&p, true, &mut no_warn).is_none());
        let bioseq = Bioseq {
            local_id: String::new(),
            title: Vec::new(),
            is_aa: true,
            length: 1,
            seq_data: None,
            parsed: Some(p),
        };
        let set = extract_bioseq_deflines(&bioseq, &mut no_warn);
        assert_eq!(
            titles_and_ids(&set),
            [(b"leading spaces".to_vec(), vec![refseq("NP_100008", 1)])]
        );
    }

    // survey B fixture p04_molmismatch and the probe of a nucleotide range.
    #[test]
    fn ids_of_the_other_molecule_become_local() {
        let mut no_warn = |_: &[u8]| panic!("no warning");
        let p = parsed(b">NC_000913.3 raw nucleotide", vec![], None);
        let set = second_parse(&p, true, &mut no_warn).unwrap();
        assert_eq!(
            set.deflines[0].ids,
            [ParsedSeqId::Local(ObjectId::Str(b"NC_000913.3".to_vec()))]
        );
        let mut ids = vec![refseq("NC_000913", 3), ParsedSeqId::Gi(5)];
        localize_other_molecule(&mut ids, true, &mut no_warn);
        assert_eq!(
            ids[0],
            ParsedSeqId::Local(ObjectId::Str(b"NC_000913.3".to_vec()))
        );
        assert_eq!(ids[1], ParsedSeqId::Gi(5));
        // A nucleotide id with a range: the Bioseq has lcl|3 and the replaced ref id in
        // its history, so the local second parse is not used.
        let mut p = parsed(
            b">NC_000913.3:1-5 t",
            vec![ParsedSeqId::Local(ObjectId::Id(3))],
            Some(b"t"),
        );
        p.hist_old_id = Some(refseq("NC_000913", 3));
        let bioseq = Bioseq {
            local_id: String::new(),
            title: Vec::new(),
            is_aa: false,
            length: 1,
            seq_data: None,
            parsed: Some(p),
        };
        let set = extract_bioseq_deflines(&bioseq, &mut no_warn);
        assert_eq!(
            titles_and_ids(&set),
            [(b"t".to_vec(), vec![ParsedSeqId::Local(ObjectId::Id(3))])]
        );
    }

    #[test]
    fn large_gis_skip_the_sequence_and_log_the_last_gi() {
        let mut no_warn = |_: &[u8]| panic!("no warning");
        assert_eq!(
            large_gi_to_skip(&[ParsedSeqId::Gi(4_294_967_295)], &mut no_warn),
            None
        );
        assert_eq!(
            large_gi_to_skip(
                &[ParsedSeqId::Gi(4_294_967_296), ParsedSeqId::Gi(7)],
                &mut no_warn
            ),
            Some(7)
        );
    }
}
