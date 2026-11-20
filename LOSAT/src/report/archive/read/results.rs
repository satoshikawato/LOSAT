//! The results of an archive object (`Blast4-get-search-results-reply`): what NCBI's archive
//! reader makes of them for blast_formatter (`CRemoteBlast::GetResultSet`): the alignments of each
//! query (`GetSeqAlignSets`), the masks of each query (`GetMasks`), the Karlin-Altschul
//! parameters (`GetKABlocks`) and the effective search space, length adjustment and matrix of the
//! search statistics (`GetSearchStats`). The reply stays the value read: the Seq-aligns are
//! written again as read by the formats that write them (`asn_types::write::write_value`), and
//! the members NCBI does not read here (`phi-alignments`, `pssm`, `simple-results`) stay in it.

use std::ops::Range;

use super::request::Queries;
use super::{seq_align_id, seq_loc_id, Stop};
use crate::ncbi_str::string_to_int8;
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::Value;
use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::ParsedSeqId;

/// The results of a search, as read.
#[derive(Debug)]
pub(crate) struct Results {
    reply: Value,
}

/// A masked range of a query (`CSeqLocInfo` of `GetMasks`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct MaskRange<'a> {
    pub(crate) id: &'a ParsedSeqId,
    pub(crate) from: i64,
    pub(crate) to: i64,
    /// `CSeqLocInfo::ETranslationFrame`: 1, 2, 3, -1, -2, -3, or 0 (`eFrameNotSet`) for a
    /// program that does not translate the query.
    pub(crate) frame: i32,
}

/// The Karlin-Altschul parameters of the results: (lambda, K, H), ungapped and gapped, 0.0 where
/// the reply has no block of the kind.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct KarlinBlocks {
    pub(crate) ungapped: (f64, f64, f64),
    pub(crate) gapped: (f64, f64, f64),
}

/// What `GetResultSet` reads of the search statistics.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct SearchStats {
    pub(crate) effective_search_space: i64,
    pub(crate) length_adjustment: i64,
    /// The matrix of an `rpsblast` or `rpstblastn` service ("Matrix: X").
    pub(crate) matrix: Option<Vec<u8>>,
}

impl Results {
    pub(crate) fn from_value(reply: Value) -> Results {
        Results { reply }
    }

    /// The reply.
    pub(crate) fn reply(&self) -> V<'_> {
        V::new(named::BLAST4_GET_SEARCH_RESULTS_REPLY, &self.reply)
    }

    /// The Seq-aligns of `alignments` (`GetAlignments()`), in order.
    pub(crate) fn alignments(&self) -> Vec<V<'_>> {
        match self.reply().get("alignments") {
            Some(set) => set.items().collect(),
            None => Vec::new(),
        }
    }

    /// `GetSeqAlignSets()`: for each query of `queries` in order (its id: `Queries::query_ids`),
    /// the alignments of that query, as a range of [`Results::alignments`]. The alignments are
    /// cut where the id of their first row changes (`CSeq_id::Compare` other than `e_YES`); each
    /// query takes the next run when its id is the run's, else no alignment. (`GetResultSet` then
    /// gives one empty set when there is no query: `alignments.resize(1)`.)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:405-468
    // ```c++
    //     	   ITERATE(CSeq_align_set::Tdata, it, al->Get()) {
    //        	   // index 0 = query, index 1 = subject
    //        	   const int query_index = 0;
    //        	   CConstRef<CSeq_id> this_id( & (*it)->GetSeq_id(query_index) );
    //
    //       	    if (current_id.Empty() || (CSeq_id::e_YES != this_id->Compare(*current_id))) {
    //               if (cur_set.NotEmpty()) {
    //                   temp.push_back(cur_set);
    //               }
    //               cur_set.Reset(new CSeq_align_set);
    //               current_id = this_id;
    //             }
    //       	    cur_set->Set().push_back(*it);
    //           }
    //     ...
    //     ITERATE(CSearchResultSet::TQueryIdVector, it,  query_ids) {
    //         const int query_index = 0;
    //         if (sap_index < temp.size())
    //         {
    //              list< CRef< CSeq_align > > sal = temp[sap_index]->Get();
    //              CConstRef<CSeq_id> this_id( & (sal.front()->GetSeq_id(query_index) ));
    //              if (CSeq_id::e_YES == (*it)->Compare(sal.front()->GetSeq_id(query_index) ))
    //              {
    //                   rv.push_back(temp[sap_index]);
    //                   sap_index++;
    //              }
    //              else
    //              {
    //                   cur_set.Reset(new CSeq_align_set);
    //                   rv.push_back(cur_set);
    //              }
    //         }
    //         else
    //         {
    //              cur_set.Reset(new CSeq_align_set);
    //              rv.push_back(cur_set);
    //         }
    //     }
    // ```
    // (A null query id, `(*it)->Compare`, is dereferenced only while a run is left.)
    pub(crate) fn alignment_sets(&self, queries: &Queries) -> Result<Vec<Range<usize>>, Stop> {
        let alignments = self.alignments();
        let mut runs: Vec<(Range<usize>, &ParsedSeqId)> = Vec::new();
        for (at, align) in alignments.iter().enumerate() {
            let id = seq_align_id(*align, 0)?;
            match runs.last_mut() {
                Some((run, current)) if compare(id, current) == Sic::Yes => run.end = at + 1,
                _ => runs.push((at..at + 1, id)),
            }
        }
        let query_ids = queries.query_ids()?;
        let mut sets = Vec::with_capacity(query_ids.len());
        let mut next = 0;
        for query in query_ids {
            match runs.get(next) {
                Some((run, first)) => {
                    let query =
                        query.ok_or(Stop::NullPointer("GetSeqAlignSets: a null query id"))?;
                    if compare(query, first) == Sic::Yes {
                        sets.push(run.clone());
                        next += 1;
                    } else {
                        sets.push(run.start..run.start);
                    }
                }
                None => sets.push(alignments.len()..alignments.len()),
            }
        }
        Ok(sets)
    }

    /// `GetMasks()`: the masked ranges of each of the `num_queries` queries
    /// (`Queries::num_queries`). The masks are taken in order: a mask whose first location's id
    /// is not the previous mask's (`CSeq_id::Compare` other than `e_YES`) belongs to the next
    /// query; each interval of that location (a `packed-int`) is a range, its frame the mask's
    /// for a program that translates the query (`query_translated`, asked only when there are
    /// masks), else none.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:489-546
    // ```c++
    //     TSeqLocInfoVector retval;
    //     retval.resize(GetQueries()->GetNumQueries());
    //
    //     TGSRR::TMasks network_masks = x_GetMasks();
    //     if (network_masks.empty()) {
    //         return retval;
    //     }
    //
    //     EBlastProgramType program = NetworkProgram2BlastProgramType(m_Program,
    //                                                                 m_Service);
    //     CConstRef<CSeq_id> previous_seqid;
    //     size_t query_index = 0;
    //
    //     ITERATE(TGSRR::TMasks, masks_for_frame, network_masks) {
    //
    //         _ASSERT(masks_for_frame->NotEmpty());
    //
    //         CConstRef<CSeq_id> current_seqid
    //             ((*masks_for_frame)->GetLocations().front()->GetId());
    //         if (previous_seqid.Empty()) {
    //             previous_seqid = current_seqid;
    //         }
    //
    //         // determine which query are we setting the masks for...
    //         TMaskedQueryRegions* mqr = NULL;
    //         if (CSeq_id::e_YES == current_seqid->Compare(*previous_seqid)) {
    //             mqr = &retval[query_index];
    //         } else {
    //             mqr = &retval[++query_index];
    //             previous_seqid = current_seqid;
    //         }
    //         ...
    //         CRef<CSeq_loc> masks =
    //             (*masks_for_frame)->GetLocations().front();
    //         _ASSERT(masks->IsPacked_int());
    //
    //         const CPacked_seqint& packed_int = masks->GetPacked_int();
    //         const EBlast4_frame_type frame = (*masks_for_frame)->GetFrame();
    //         ITERATE(CPacked_seqint::Tdata, mask, packed_int.Get()) {
    //             CRef<CSeq_interval> si
    //                 (new CSeq_interval(const_cast<CSeq_id&>((*mask)->GetId()),
    //                                    (*mask)->GetFrom(), (*mask)->GetTo()));
    //             CRef<CSeqLocInfo> sli
    //                 (new CSeqLocInfo(si, NetworkFrame2FrameNumber(frame, program)));
    //             mqr->push_back(sli);
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1990-2017
    // ```c++
    // CSeqLocInfo::ETranslationFrame
    // NetworkFrame2FrameNumber(objects::EBlast4_frame_type frame,
    //                          EBlastProgramType program)
    // {
    //     if (Blast_QueryIsTranslated(program)) {
    //         switch (frame) {
    //         case eBlast4_frame_type_plus1:  return CSeqLocInfo::eFramePlus1;
    //         case eBlast4_frame_type_plus2:  return CSeqLocInfo::eFramePlus2;
    //         case eBlast4_frame_type_plus3:  return CSeqLocInfo::eFramePlus3;
    //         case eBlast4_frame_type_minus1: return CSeqLocInfo::eFrameMinus1;
    //         case eBlast4_frame_type_minus2: return CSeqLocInfo::eFrameMinus2;
    //         case eBlast4_frame_type_minus3: return CSeqLocInfo::eFrameMinus3;
    //         default: abort();
    //         }
    //         _TROUBLE;
    //     }
    //     ...
    //     return CSeqLocInfo::eFrameNotSet;
    // }
    // ```
    // (`_ASSERT`s are not in NCBI's release build. A mask without locations (`front()` of an
    // empty list), or one more query than `GetNumQueries` (`retval[++query_index]`), reads
    // memory the input does not decide; a first location without one id is a null `CConstRef`.)
    pub(crate) fn masks(
        &self,
        num_queries: usize,
        query_translated: impl FnOnce() -> Result<bool, Stop>,
    ) -> Result<Vec<Vec<MaskRange<'_>>>, Stop> {
        let mut out: Vec<Vec<MaskRange<'_>>> = vec![Vec::new(); num_queries];
        let Some(masks) = self
            .reply()
            .get("masks")
            .filter(|masks| masks.item(0).is_some())
        else {
            return Ok(out);
        };
        let translated = query_translated()?;
        let mut previous: Option<&ParsedSeqId> = None;
        let mut index = 0;
        for mask in masks.items() {
            let first = mask.field("locations").item(0).ok_or(Stop::Undefined(
                "GetMasks: front() of a mask without locations",
            ))?;
            let current = seq_loc_id(first);
            if previous.is_none() {
                previous = current;
            }
            let null = Stop::NullPointer("GetMasks: the id of a mask's first location");
            let (current, last) = (current.ok_or(null.clone())?, previous.ok_or(null)?);
            if compare(current, last) != Sic::Yes {
                index += 1;
                previous = Some(current);
            }
            let ranges = out.get_mut(index).ok_or(Stop::Undefined(
                "GetMasks: a mask of a query after the last (retval[++query_index])",
            ))?;
            let (variant, intervals) = first.choice();
            if variant != "packed-int" {
                return Err(Stop::InvalidSelection {
                    module: "NCBI-Seqloc",
                    type_name: "Seq-loc",
                    required: "packed-int",
                    current: variant,
                });
            }
            let frame = mask.field("frame").int();
            for interval in intervals.items() {
                let frame = match (translated, frame) {
                    (false, _) => 0,
                    (true, 1..=3) => frame as i32,
                    (true, 4..=6) => 3 - frame as i32,
                    (true, _) => {
                        return Err(Stop::Abort(
                            "NetworkFrame2FrameNumber: a frame other than plus1-minus3",
                        ))
                    }
                };
                ranges.push(MaskRange {
                    id: interval.field("id").seq_id(),
                    from: interval.field("from").int(),
                    to: interval.field("to").int(),
                    frame,
                });
            }
        }
        Ok(out)
    }

    /// The Karlin-Altschul blocks as `GetResultSet` takes them: in order, a gapped block sets
    /// the gapped parameters and another the ungapped ones, until both kinds are found.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2218-2238
    // ```c++
    //         bool found_gapped = false, found_ungapped = false;
    //         pair<double, double> lambdas, Ks, Hs;
    //         TKarlinAltschulBlocks ka_blocks = GetKABlocks();
    //
    //         ITERATE(TKarlinAltschulBlocks, itr, ka_blocks) {
    //             if ((*itr)->GetGapped()) {
    //                 lambdas.second = (*itr)->GetLambda();
    //                 Ks.second = (*itr)->GetK();
    //                 Hs.second = (*itr)->GetH();
    //                 found_gapped = true;
    //             } else {
    //                 lambdas.first = (*itr)->GetLambda();
    //                 Ks.first = (*itr)->GetK();
    //                 Hs.first = (*itr)->GetH();
    //                 found_ungapped = true;
    //             }
    //
    //             if (found_gapped && found_ungapped) {
    //                 break;
    //             }
    //         }
    // ```
    pub(crate) fn karlin_blocks(&self) -> KarlinBlocks {
        let mut out = KarlinBlocks::default();
        let (mut gapped, mut ungapped) = (false, false);
        let Some(blocks) = self.reply().get("ka-blocks") else {
            return out;
        };
        for block in blocks.items() {
            let real = |name: &str| match block.field(name).value {
                Value::Real(real) => *real,
                other => unreachable!("a REAL of the table, not {other:?}"),
            };
            let params = (real("lambda"), real("k"), real("h"));
            if matches!(block.field("gapped").value, Value::Bool(true)) {
                out.gapped = params;
                gapped = true;
            } else {
                out.ungapped = params;
                ungapped = true;
            }
            if gapped && ungapped {
                break;
            }
        }
        out
    }

    /// The search statistics as `GetResultSet` reads them for the `service` of the request: in
    /// each string, "Effective search space used: " and then "Length adjustment: " are removed
    /// where found and the rest read as a number (0 when it is not one; the last string that
    /// has one gives it); for `rpsblast` and `rpstblastn`, the first string (as left) with
    /// "Matrix: " gives the matrix, the rest after the text is removed.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2187-2215
    // ```c++
    //         const string kTarget("Effective search space used: ");
    //         const string kLength_Adj("Length adjustment: ");
    //         list<string> search_stats = GetSearchStats();
    //         Int8 effective_search_space = 0;
    //         Int8 length_adjustment = 0;
    //         NON_CONST_ITERATE(list<string>, itr, search_stats) {
    //             if (NStr::Find(*itr, kTarget) != NPOS) {
    //                 NStr::ReplaceInPlace(*itr, kTarget, kEmptyStr);
    //                 effective_search_space =
    //                     NStr::StringToInt8(*itr, NStr::fConvErr_NoThrow);
    //             }
    //             if (NStr::Find(*itr, kLength_Adj) != NPOS) {
    //                             NStr::ReplaceInPlace(*itr, kLength_Adj, kEmptyStr);
    //                             length_adjustment =
    //                                 NStr::StringToInt8(*itr, NStr::fConvErr_NoThrow);
    //             }
    //         }
    //
    //         if((m_Service == "rpsblast") || (m_Service == "rpstblastn") )
    //         {
    //         	const string kTmp("Matrix: ");
    //         	 NON_CONST_ITERATE(list<string>, itr, search_stats) {
    //         		 if (NStr::Find(*itr, kTmp) != NPOS) {
    //         			 NStr::ReplaceInPlace(*itr, kTmp, kEmptyStr);
    //         			 m_CBOH->SetOptions().SetMatrixName((*itr).c_str());
    //         			 break;
    //         		 }
    //         	 }
    //         }
    // ```
    pub(crate) fn search_stats(&self, service: &[u8]) -> SearchStats {
        let mut out = SearchStats::default();
        let mut stats: Vec<Vec<u8>> = match self.reply().get("search-stats") {
            Some(list) => list.items().map(|text| text.str().to_vec()).collect(),
            None => return out,
        };
        let number = |text: &[u8]| string_to_int8(text, 0).unwrap_or(0);
        for text in &mut stats {
            if replace_all(text, b"Effective search space used: ") {
                out.effective_search_space = number(text);
            }
            if replace_all(text, b"Length adjustment: ") {
                out.length_adjustment = number(text);
            }
        }
        if service == b"rpsblast" || service == b"rpstblastn" {
            for mut text in stats {
                if replace_all(&mut text, b"Matrix: ") {
                    out.matrix = Some(text);
                    break;
                }
            }
        }
        out
    }
}

/// `NStr::Find(text, search) != NPOS`, then `NStr::ReplaceInPlace(text, search, kEmptyStr)`:
/// every `search` removed, each search going on from where the last was removed. Whether one
/// was found.
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3401-3428
// ```c++
// string& NStr::ReplaceInPlace(string& src,
//                              const string& search, const string& replace,
//                              SIZE_TYPE start_pos, SIZE_TYPE max_replace,
//                              SIZE_TYPE* num_replace)
// {
//     ...
//     if ( start_pos + search.size() > src.size()  ||  search == replace )
//         return src;
//
//     bool equal_len = (search.size() == replace.size());
//     for (SIZE_TYPE count = 0; !(max_replace && count >= max_replace); count++){
//         start_pos = src.find(search, start_pos);
//         if (start_pos == NPOS)
//             break;
//         ...
//             src.replace(start_pos, search.size(), replace);
//         ...
//         start_pos += replace.size();
//         ...
//     }
//     return src;
// }
// ```
// (`start_pos` 0 and `max_replace` 0 by default; the search was found, so the first test
// does not return.)
fn replace_all(text: &mut Vec<u8>, search: &[u8]) -> bool {
    let mut from = 0;
    let mut found = false;
    while let Some(at) = text[from..]
        .windows(search.len())
        .position(|window| window == search)
    {
        let at = from + at;
        text.drain(at..at + search.len());
        from = at;
        found = true;
    }
    found
}
