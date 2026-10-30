//! `-outfmt %f` (the default): NCBI's `CBlastDB_FastaFormatter`. The Bioseq of a sequence
//! (`CSeqDB::GetBioseq`: the ids of its first definition line, a title joined from every
//! definition line), the header of `CBlastDeflineUtil::ProcessFastaDeflines` (short ids) or
//! of `CFastaOstream` (`-long_seqids`), and the sequence lines of
//! `CFastaOstream::WriteSequence`.

use std::collections::HashSet;
use std::io::Write;

use crate::report::defline::{bioseq_title_from_source, generate_bioseq_defline};
use crate::seqid::accguide::parse_flags::PARTIAL_OK;
use crate::seqid::label::{bare, fasta};
use crate::seqid::model::{Dbtag, ObjectId};
use crate::seqid::parse::parse_ids;
use crate::seqid::rank::{fasta_aa_score, fasta_na_score, find_best, text_score};
use crate::seqid::{DeflineSet, ParsedSeqId};

use super::letters::reverse_complement;
use super::retrieve::{numeric_target, parse_target, Config, Source, Stop};

/// `CBlastDB_FastaFormatter` with its one `CFastaOstream` (`m_fasta`), whose record of the
/// ids written keeps `-long_seqids` from writing one twice.
pub(super) struct FastaFormatter {
    width: u32,
    long_seqids: bool,
    /// `CFastaOstream::m_PreviousWholeIds`: the FASTA strings of the ids written without a
    /// range.
    previous_whole_ids: HashSet<Vec<u8>>,
}

impl FastaFormatter {
    pub fn new(width: u32, long_seqids: bool) -> Self {
        FastaFormatter {
            width,
            long_seqids,
            previous_whole_ids: HashSet::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:436-531
    // ```c++
    // int CBlastDB_FastaFormatter::Write(CSeqDB::TOID oid, const CBlastDB_FormatterConfig & config, string target_id)
    // {
    // 	int status = -1;
    // 	CRef<CBioseq> bioseq = m_BlastDb.GetBioseq(oid);
    // 	if(target_id != kEmptyStr) {
    // 		Int8 num_id = NStr::StringToNumeric<Int8>(target_id, NStr::fConvErr_NoThrow);
    // 		if(errno) {
    // 			CSeq_id seq_id(target_id, CSeq_id::fParse_PartialOK | CSeq_id::fParse_Default);
    // 			bioseq = m_BlastDb.GetBioseq(oid, ZERO_GI, &seq_id);
    // 		}
    // 		else {
    // 			bioseq = m_BlastDb.GetBioseq(oid, GI_FROM(Int8, num_id));
    // 		}
    // 	}
    // 	...
    // 	if (bioseq.Empty()) {
    // 		return status;
    // 	}
    // 	...
    //     CRef<CSeq_loc> range;
    //     if(config.m_SeqRange.NotEmpty()) {
    //     	TSeqPos length = m_BlastDb.GetSeqLength(oid);
    //     	TSeqRange r = config.m_SeqRange;
    //         if((TSeqPos)length <= config.m_SeqRange.GetTo()) {
    //         	r.SetTo(length-1);
    //     	}
    //         if(!r.Empty()) {
    //         	range.Reset(new CSeq_loc(*(FindBestChoice(bioseq->GetId(), CSeq_id::BestRank)), r.GetFrom(), r.GetTo()));
    //         }
    //     }
    //     if(m_UseLongSeqIds) {
    //     	...
    //     }
    //     else {
    //     	string fasta_deflines = kEmptyStr;
    //     	...
    //     	m_Out << fasta_deflines;
    //     	...
    //     	m_fasta.WriteSequence(scope.AddBioseq(*bioseq), range.GetPointer());
    //    }
    //    return 0;
    // }
    // ```
    /// Writes the FASTA of sequence `oid` to `out` (nothing for a sequence NCBI's `GetBioseq`
    /// gives no Bioseq for: no definition line, or no letters).
    pub fn write(
        &mut self,
        source: &Source,
        oid: u32,
        config: &Config,
        target: Option<&str>,
        out: &mut dyn Write,
    ) -> Result<(), Stop> {
        let set = source.deflines(oid)?;
        let set = match target {
            Some(target) => target_set(set, target)?,
            None => set,
        };
        let Some(defline) = set.deflines.first() else {
            return Ok(());
        };
        if defline.ids.is_empty() {
            return Err(Stop::Refused(format!(
                "a definition line without Seq-ids is not supported by LOSAT's {}",
                super::PROGRAM
            )));
        }
        let length = source.length(oid)?;
        if length < 1 {
            return Ok(());
        }
        let range = config.range.clipped(length);
        let mut text = Vec::new();
        if self.long_seqids {
            self.long_header(source, &set, range, config, &mut text)?;
        } else {
            short_header(&set, range, config, &mut text)?;
        }
        let mut letters = source.fasta_letters(oid, range)?;
        if config.strand_minus {
            reverse_complement(&mut letters);
        }
        write_sequence(&letters, self.width, &mut text);
        out.write_all(&text).map_err(Stop::io)
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2927-2975
    // ```c++
    // void CFastaOstream::x_WriteSeqIds(const CBioseq& bioseq,
    //                                   const CSeq_loc* location)
    // {
    //     bool have_range = (location != NULL  &&  !location->IsWhole()
    //                        &&  !(m_Flags & fSuppressRange) );
    //
    //     if ( !have_range  &&  (m_Flags & fNoDupCheck) == 0) {
    //         ITERATE (CBioseq::TId, id, bioseq.GetId()) {
    //             CSeq_id_Handle idh = CSeq_id_Handle::GetHandle(**id);
    //             pair<TSeq_id_HandleSet::iterator, bool> p
    //                 = m_PreviousWholeIds.insert(idh);
    //             if ( !p.second ) {
    //                 NCBI_THROW(CObjmgrUtilException, eBadLocation,
    //                            "Duplicate Seq-id " + (*id)->AsFastaString()
    //                            + " in FASTA output");
    //             }
    //         }
    //     }
    //
    //     m_Out << '>';
    //     ...
    //         x_WriteAsFasta(bioseq);
    //     ...
    //     if (have_range) {
    //         char delim = ':';
    //         for (CSeq_loc_CI it(*location);  it;  ++it) {
    //             CSeq_loc::TRange range = it.GetRange();
    //             TSeqPos from = range.GetFrom() + 1, to = range.GetTo() + 1;
    //             _ASSERT(from <= to);
    //             m_Out << delim;
    //             if (it.IsSetStrand()  &&  IsReverse(it.GetStrand())) {
    //                 m_Out << 'c' << to << '-' << from;
    //             } else {
    //                 m_Out << from << '-' << to;
    //             }
    //             delim = ',';
    //         }
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2875-2925
    // ```c++
    //     x_GetBestId(gi_id, best_id, hide_prefix, bioseq);
    //
    //     if (best_id.NotEmpty())
    //     {
    //         // RW-139, no GI in FASTA output
    //         if (gi_id.NotEmpty() && (m_Flags & fEnableGI) && !best_id->IsGi())
    //         {
    //             gi_id->WriteAsFasta(m_Out);
    //             m_Out << '|';
    //         }
    //         ...
    //             best_id->WriteAsFasta(m_Out);
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3000-3020
    // ```c++
    //     string safe_title = (!custom_title.empty()) ? custom_title
    //         : m_Gen->GenerateDefline(bioseq_handle, x_GetTitleFlags());
    //
    //     if ( !safe_title.empty() ) {
    //         ...
    //         if (safe_title[0] != ' ') {
    //             m_Out << ' ';
    //         }
    //         ...
    //         m_Out << safe_title;
    //     }
    //     m_Out << '\n';
    // ```
    /// The header of `-long_seqids`: the Bioseq's best id for FASTA (`FastaNARank`,
    /// `FastaAARank`) after its GI, the range as `:a-b` (`:cb-a` for the minus strand), and
    /// the title `CDeflineGenerator` makes of the joined title (line BDC DECISIONS 12: NCBI
    /// writes the forward letters under the `c` header of a range on the minus strand;
    /// LOSAT writes the reverse complement, as the short header's path does).
    fn long_header(
        &mut self,
        source: &Source,
        set: &DeflineSet,
        range: Option<(u32, u32)>,
        config: &Config,
        out: &mut Vec<u8>,
    ) -> Result<(), Stop> {
        let ids = &set.deflines[0].ids;
        if range.is_none() {
            for id in ids {
                let mut key = Vec::new();
                fasta(id, &mut key);
                if !self.previous_whole_ids.insert(key.clone()) {
                    return Err(Stop::Ncbi(format!(
                        "Duplicate Seq-id {} in FASTA output",
                        String::from_utf8_lossy(&key)
                    )));
                }
            }
        }
        out.push(b'>');
        let best = if source.protein {
            find_best(ids, fasta_aa_score)
        } else {
            find_best(ids, fasta_na_score)
        }
        .expect("the definition line has ids");
        if let Some(gi) = ids.iter().find(|id| id.gi().is_some()) {
            if best.gi().is_none() {
                fasta(gi, out);
                out.push(b'|');
            }
        }
        fasta(best, out);
        if let Some((from, to)) = range {
            if config.strand_minus {
                out.extend_from_slice(format!(":c{}-{}", to + 1, from + 1).as_bytes());
            } else {
                out.extend_from_slice(format!(":{}-{}", from + 1, to + 1).as_bytes());
            }
        }
        let mut title = bioseq_title(set);
        if config.ctrl_a {
            title = rejoin_defline_title(&title, true)?;
        }
        let title = generated_title(source, set, &title)?;
        if !title.is_empty() {
            if title[0] != b' ' {
                out.push(b' ');
            }
            out.extend_from_slice(&title);
        }
        out.push(b'\n');
        Ok(())
    }
}

/// The definition lines NCBI's `GetBioseq` keeps for `-target_only`: the first that holds
/// the target id (`s_SeqDB_SeqIdIn`), or the local id of its label.
fn target_set(set: DeflineSet, target: &str) -> Result<DeflineSet, Stop> {
    let target_id = match numeric_target(target) {
        Some(number) => ParsedSeqId::Gi(number),
        None => parse_target(target)?,
    };
    let Some(index) = set.bioseq_defline(&target_id) else {
        return Err(Stop::Ncbi(
            "Error: oid headers do not contain target gi/seq_id.".to_string(),
        ));
    };
    Ok(DeflineSet {
        deflines: vec![set.deflines[index].clone()],
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:957-1019
// ```c++
// static void
// s_GetBioseqTitle(CRef<CBlast_def_line_set> deflines, string & title)
// {
//     title.erase();
//     ...
//     bool first_defline(true);
//
//     for(TDefIt iter = dl.begin(); iter != dl.end(); iter++) {
//         ostringstream oss;
//
//         const CBlast_def_line & defline = **iter;
//
//         if (! title.empty()) {
//             //oss << "\1";
//             oss << " ";
//         }
//
//         bool wrote_seqids(false);
//
//         if ((!first_defline) && defline.CanGetSeqid()) {
//             const list< CRef<CSeq_id > > & sl = defline.GetSeqid();
//
//             bool first_seqid(true);
//
//             // First should look like: "<title>"
//             // Others should look like: " ><seqid>|<seqid>|<seqid><title>"
//
//             for (TSeqIt seqit = sl.begin(); seqit != sl.end(); seqit++) {
//                 if (first_seqid) {
//                     oss << ">";
//                 } else {
//                     oss << "|";
//                 }
//
//                 (*seqit)->WriteAsFasta(oss);
//
//                 first_seqid = false;
//                 wrote_seqids = true;
//             }
//         }
//
//         // Omit seqids from first defline
//         first_defline = false;
//
//         if (defline.CanGetTitle()) {
//             if (wrote_seqids) {
//                 oss << " ";
//             }
//             oss << defline.GetTitle();
//         }
//
//         title += oss.str();
//     }
// }
// ```
/// The title of the Bioseq: the first definition line's title, then ` >` with the FASTA
/// strings of the ids of each later one (joined by `|`) and its title.
pub(super) fn bioseq_title(set: &DeflineSet) -> Vec<u8> {
    let mut title = Vec::new();
    let mut first_defline = true;
    for defline in &set.deflines {
        let mut oss = Vec::new();
        if !title.is_empty() {
            oss.push(b' ');
        }
        let mut wrote_seqids = false;
        if !first_defline {
            for (at, id) in defline.ids.iter().enumerate() {
                oss.push(if at == 0 { b'>' } else { b'|' });
                fasta(id, &mut oss);
                wrote_seqids = true;
            }
        }
        first_defline = false;
        if let Some(text) = &defline.title {
            if wrote_seqids {
                oss.push(b' ');
            }
            oss.extend_from_slice(text);
        }
        title.extend_from_slice(&oss);
    }
    title
}

/// The ids `CSeq_id::ParseIDs(ids, text, fParse_PartialOK)` gives (none when it throws).
/// `Err`: a parse that posts a warning NCBI prints (an id of several parts of which one
/// fails), whose line LOSAT does not reproduce; blastdbcmd filters out the warning of extra
/// parts (`CSeq_id` 1306.10) alone.
fn parse_ids_partial(text: &[u8]) -> Result<Vec<ParsedSeqId>, Stop> {
    let mut warned = None;
    let parsed = parse_ids(text, PARTIAL_OK, &mut |warning| {
        if warning.subcode != 10 {
            warned = Some(warning.message);
        }
    });
    if let Some(message) = warned {
        return Err(Stop::Refused(format!(
            "a Seq-id in a title that NCBI reads with a warning ({message}) is not supported \
             by LOSAT's {}",
            super::PROGRAM
        )));
    }
    Ok(parsed
        .map(|ids| ids.into_iter().flatten().collect())
        .unwrap_or_default())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:630-665
// ```c++
// static string
// s_ConfigureDeflineTitle(const string& title, bool use_ctrl_a)
// {
//     static const string kStandardSeparator(" >");
//     const string kSeparator(use_ctrl_a ? CTRL_A : kStandardSeparator);
//     string retval;
//     list<string> tokens;
//     NStr::Split(title, kStandardSeparator, tokens, NStr::fSplit_ByPattern);
//     int idx = 0;
//     for (auto token : tokens) {
//         if (idx++ == 0) {
//             retval += token;
//             continue;
//         }
//         SIZE_TYPE pos = token.find(' ');
//         const string kPossibleId(token, 0, pos != NPOS ? pos : token.length());
//         CBioseq::TId seqids;
//
//         try {
//             CSeq_id::ParseIDs(seqids, kPossibleId, CSeq_id::fParse_PartialOK);
//         } catch (const CException&) {}
//
//         if (!seqids.empty()) {
//             retval += kSeparator;
//             CRef<CSeq_id> id = FindBestChoice(seqids, CSeq_id::Score);
//             retval += GetBareId(*id);
//             if (pos != NPOS)
//                 retval += token.substr(pos, token.length() - pos);
//         } else {
//             retval += kStandardSeparator + token;
//         }
//     }
//     return retval;
// }
// ```
/// The title of the short header: each ` >` part whose first word parses as Seq-ids gets
/// the bare form of their best, after `\x01` with `-ctrl_a`.
fn configure_defline_title(title: &[u8], use_ctrl_a: bool) -> Result<Vec<u8>, Stop> {
    let separator: &[u8] = if use_ctrl_a { b"\x01" } else { b" >" };
    let mut retval = Vec::new();
    for (idx, token) in split_by_pattern(title, b" >").into_iter().enumerate() {
        if idx == 0 {
            retval.extend_from_slice(token);
            continue;
        }
        let pos = token.iter().position(|&c| c == b' ');
        let possible_id = &token[..pos.unwrap_or(token.len())];
        let seqids = parse_ids_partial(possible_id)?;
        match find_best(&seqids, text_score) {
            Some(id) => {
                retval.extend_from_slice(separator);
                bare(id, &mut retval);
                if let Some(pos) = pos {
                    retval.extend_from_slice(&token[pos..]);
                }
            }
            None => {
                retval.extend_from_slice(b" >");
                retval.extend_from_slice(token);
            }
        }
    }
    Ok(retval)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:362-398
// ```c++
// static string
// s_RejoinDeflineTitle(const string& title, bool use_ctrl_a )
// {
//     static const string kStandardSeparator(" >");
//     const string kSeparator(use_ctrl_a ? "\001" : kStandardSeparator); // actual, will b eused for output
//     string retval;
//     list<string> tokens;
//     NStr::Split(title, kStandardSeparator, tokens, NStr::fSplit_ByPattern);
//     int idx = 0;
//     for (auto token : tokens) {
//         if (idx++ == 0) {
//             retval += token;
//             continue;
//         }
//
//         SIZE_TYPE pos = token.find(' ');
//         const string kPossibleId(token, 0, pos != NPOS ? pos : token.length());
//         CBioseq::TId seqids;
// 	// verify found token via CSeq_id::ParseIDs
//         try {
//             CSeq_id::ParseIDs(seqids, kPossibleId, CSeq_id::fParse_PartialOK);
//         } catch (const CException&) {}
//
//         if (!seqids.empty()) {
//             retval += kSeparator + token;
//         } else {
//             retval += kStandardSeparator + token;
//         }
//
//     }
//     return retval;
// }
// ```
/// The title of `-long_seqids -ctrl_a`: `\x01` before each ` >` part whose first word
/// parses as Seq-ids, the part itself unchanged.
fn rejoin_defline_title(title: &[u8], use_ctrl_a: bool) -> Result<Vec<u8>, Stop> {
    let separator: &[u8] = if use_ctrl_a { b"\x01" } else { b" >" };
    let mut retval = Vec::new();
    for (idx, token) in split_by_pattern(title, b" >").into_iter().enumerate() {
        if idx == 0 {
            retval.extend_from_slice(token);
            continue;
        }
        let pos = token.iter().position(|&c| c == b' ');
        let possible_id = &token[..pos.unwrap_or(token.len())];
        if parse_ids_partial(possible_id)?.is_empty() {
            retval.extend_from_slice(b" >");
        } else {
            retval.extend_from_slice(separator);
        }
        retval.extend_from_slice(token);
    }
    Ok(retval)
}

/// `NStr::Split(text, pattern, tokens, fSplit_ByPattern)`: the pieces between the
/// occurrences of the whole pattern, empty pieces kept.
fn split_by_pattern<'a>(text: &'a [u8], pattern: &[u8]) -> Vec<&'a [u8]> {
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut at = 0;
    while at + pattern.len() <= text.len() {
        if &text[at..at + pattern.len()] == pattern {
            tokens.push(&text[start..at]);
            at += pattern.len();
            start = at;
        } else {
            at += 1;
        }
    }
    tokens.push(&text[start..]);
    tokens
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:999-1032
// ```c++
//     string range;
//     if (location != NULL) {
//         TSeqPos start = location->GetStart(eExtreme_Biological) + 1;
//         TSeqPos stop = location->GetStop(eExtreme_Biological) + 1;
//         if (strand == eNa_strand_minus) {
//             range = ":c"
//                     + NStr::IntToString(stop) + "-" + NStr::IntToString(start)
//                     + " ";
//         } else {
//             range = ":"
//                     + NStr::IntToString(start) + "-" + NStr::IntToString(stop)
//                     + " ";
//         }
//     }
//     if (id->IsGeneral() && id->GetGeneral().GetDb() == "BL_ORD_ID") {
//         out = ">" + range + s_GetTitle(bioseq) + '\n';
//     }
//     else if (id->IsLocal()) {
//         string lcl_tmp = id->AsFastaString();
//         lcl_tmp = lcl_tmp.erase(0,4);
//         out = ">" + lcl_tmp + (range.empty() ? " " : range)
//                 + s_GetTitle(bioseq) + '\n';
//     } else {
//         out = '>';
//         id = FindBestChoice(bioseq.GetId(), CSeq_id::Score);
//         out += GetBareId(*id) + (range.empty() ? " " : range);
//
//         string title = s_GetTitle(bioseq);
//         out += s_ConfigureDeflineTitle(title, use_ctrla);
//         out += '\n';
//     }
// ```
/// The header of the short ids (`ProcessFastaDeflines`): the title alone after
/// `gnl|BL_ORD_ID|n`, the local id's string, or the bare form of the Bioseq's best id
/// (`Score`) and the configured title; the range as `:a-b ` or `:cb-a ` in the place of the
/// space.
fn short_header(
    set: &DeflineSet,
    range: Option<(u32, u32)>,
    config: &Config,
    out: &mut Vec<u8>,
) -> Result<(), Stop> {
    let ids = &set.deflines[0].ids;
    let first = &ids[0];
    let range_text = range.map(|(from, to)| {
        if config.strand_minus {
            format!(":c{}-{} ", to + 1, from + 1)
        } else {
            format!(":{}-{} ", from + 1, to + 1)
        }
    });
    let title = bioseq_title(set);
    out.push(b'>');
    match first {
        ParsedSeqId::General(Dbtag { db, .. }) if db == b"BL_ORD_ID" => {
            if let Some(text) = &range_text {
                out.extend_from_slice(text.as_bytes());
            }
            out.extend_from_slice(&title);
        }
        ParsedSeqId::Local(_) => {
            let mut text = Vec::new();
            fasta(first, &mut text);
            out.extend_from_slice(&text[4..]);
            out.extend_from_slice(range_text.as_deref().unwrap_or(" ").as_bytes());
            out.extend_from_slice(&title);
        }
        _ => {
            let best = find_best(ids, text_score).expect("the definition line has ids");
            bare(best, out);
            out.extend_from_slice(range_text.as_deref().unwrap_or(" ").as_bytes());
            out.extend_from_slice(&configure_defline_title(&title, config.ctrl_a)?);
        }
    }
    out.push(b'\n');
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3090-3110
// ```c++
// static size_t s_TitleEndsInOrganism (
//     string& title,
//     CTempString taxname
// )
//
// {
//     size_t  pos, idx;
//     size_t  len1, len2;
//
//     len1 = title.length();
//     len2 = taxname.length();
//
//     idx = len1 - len2 - 3;
//     if (len1 > len2 + 4 && title [idx] == ' ' && title [idx + 1] == '[' && title [len1 - 1] == ']') {
//         pos = NStr::Find(title, taxname, NStr::eNocase, NStr::eReverseSearch);
//         if (pos == idx + 2) {
//             return pos - 1;
//         }
//     }
//
//     return NPOS;
// }
// ```
fn title_ends_in_organism(title: &[u8], taxname: &[u8]) -> Option<usize> {
    let (len1, len2) = (title.len(), taxname.len());
    if len1 > len2 + 4 {
        let idx = len1 - len2 - 3;
        if title[idx] == b' ' && title[idx + 1] == b'[' && title[len1 - 1] == b']' {
            let pos = (0..=len1 - len2)
                .rev()
                .find(|&at| title[at..at + len2].eq_ignore_ascii_case(taxname));
            if pos == Some(idx + 2) {
                return Some(idx + 1);
            }
        }
    }
    None
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3049-3055
// ```c++
// static inline void s_TrimMainTitle (string& str)
// {
//     size_t pos = str.find_last_not_of (".,;~ ");
//     if (pos != NPOS) {
//         str.erase (pos + 1);
//     }
// }
// ```
fn trim_main_title(title: &mut Vec<u8>) {
    if let Some(pos) = title.iter().rposition(|c| !b".,;~ ".contains(c)) {
        title.truncate(pos + 1);
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:99-133
// ```c++
// static const char* s_proteinOrganellePrefix [] = {
//   "",
//   "",
//   "chloroplast",
//   "chromoplast",
//   "kinetoplast",
//   "mitochondrion",
//   "plastid",
//   "macronuclear",
//   "extrachromosomal",
//   "plasmid",
//   "",
//   "",
//   "cyanelle",
//   "proviral",
//   "virus",
//   "nucleomorph",
//   "apicoplast",
//   "leucoplast",
//   "protoplast",
//   "endogenous virus",
//   "hydrogenosome",
//   "chromosome",
//   "chromatophore",
// };
// ```
/// `s_proteinOrganellePrefix` from `chloroplast` to `chromatophore`, the non-empty ones.
const PROTEIN_ORGANELLE_PREFIXES: [&[u8]; 19] = [
    b"chloroplast",
    b"chromoplast",
    b"kinetoplast",
    b"mitochondrion",
    b"plastid",
    b"macronuclear",
    b"extrachromosomal",
    b"plasmid",
    b"cyanelle",
    b"proviral",
    b"virus",
    b"nucleomorph",
    b"apicoplast",
    b"leucoplast",
    b"protoplast",
    b"endogenous virus",
    b"hydrogenosome",
    b"chromosome",
    b"chromatophore",
];

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3263-3410
// ```c++
// void CDeflineGenerator::x_AdjustProteinTitleSuffix (
//     const CBioseq_Handle& bsh
// )
//
// {
//     ...
//     if (m_Source.Empty()) return;
//
//     if (m_Source->IsSetTaxname()) {
//         m_Taxname = m_Source->GetTaxname();
//     }
//     ...
//     s_TrimMainTitle (m_MainTitle);
//
//     len1 = (int) m_MainTitle.length();
//     len2 = (int) m_Taxname.length();
//
//     // find [taxname]
//
//     if (len1 > len2 + 4) {
//         tpos = s_TitleEndsInOrganism(m_MainTitle, m_Taxname);
//         ...
//     }
//
//     /* do not change unless [genus species] was at the end */
//     if (tpos == NPOS) return;
//
//     m_MainTitle.erase (tpos);
//     s_TrimMainTitle (m_MainTitle);
//     len1 = (int) m_MainTitle.length();
//
//     // find (organelle)
//
//     if (len1 > 2 && m_MainTitle [len1 - 1] == ')') {
//         pos = m_MainTitle.find_last_of ("(");
//         if (pos != NPOS) {
//             for ( genome = NCBI_GENOME(chloroplast); genome <= NCBI_GENOME(chromatophore); genome++ ) {
//                 string str = s_proteinOrganellePrefix [genome];
//                 if ( ! str.empty() ) {
//                     string paren = "(" + str + ")";
//                     if (NStr::EndsWith (m_MainTitle, paren )) {
//                         opos = pos;
//                         break;
//                     }
//                 }
//             }
//         }
//         s_TrimMainTitle (m_MainTitle);
//         len1 = (int) m_MainTitle.length();
//     }
//
//     if (opos != NPOS) {
//         m_MainTitle.erase (opos);
//         s_TrimMainTitle (m_MainTitle);
//         len1 = (int) m_MainTitle.length();
//     }
//
//     if ( NStr::EndsWith (m_MainTitle, ", partial")) {
//         m_MainTitle.erase(m_MainTitle.length() - 9);
//         s_TrimMainTitle (m_MainTitle);
//     }
//
//     // then reconstruct partial (organelle) [taxname] suffix
//
//     if ( !x_IsComplete()) {
//         m_MainTitle += ", partial";
//     }
//
//     if (m_OmitTaxonomicName) return;
//
//     CTempString taxname = m_Taxname;
//
//     if (m_Genome >= NCBI_GENOME(chloroplast) && m_Genome <= NCBI_GENOME(chromatophore)) {
//         ...
//     }
//
//     // check for special taxname, go to overlapping source feature
//     ...
//     if (m_IsCrossKingdom && ! m_FirstSuperKingdom.empty() && ! m_SecondSuperKingdom.empty()) {
//         m_MainTitle += " [" + string(m_FirstSuperKingdom) + "][" + string(m_SecondSuperKingdom) + "]";
//     } else if (! taxname.empty() /* && m_MainTitle.find(taxname) == NPOS */) {
//         m_MainTitle += " [" + string(taxname) + "]";
//     }
// }
// ```
/// `x_AdjustProteinTitleSuffix` for the Bioseq of a database protein whose source has the
/// scientific name `taxname` (the taxonomy database's name of the first definition line
/// with a taxid; `GetBioseq` gives the Org-ref no name, genome or binomial without it): a
/// title that ends in ` [taxname]` loses that, a `(organelle)` and `, partial` before it,
/// and gets ` [taxname]` back (the Bioseq has no MolInfo: it is complete).
fn adjust_protein_title_suffix(title: &mut Vec<u8>, taxname: &[u8]) {
    trim_main_title(title);
    let Some(tpos) = title_ends_in_organism(title, taxname) else {
        return;
    };
    title.truncate(tpos);
    trim_main_title(title);
    let mut opos = None;
    if title.len() > 2 && title[title.len() - 1] == b')' {
        if let Some(pos) = title.iter().rposition(|&c| c == b'(') {
            for prefix in PROTEIN_ORGANELLE_PREFIXES {
                let paren = [b"(", prefix, b")"].concat();
                if title.ends_with(&paren) {
                    opos = Some(pos);
                    break;
                }
            }
        }
        trim_main_title(title);
    }
    if let Some(opos) = opos {
        title.truncate(opos);
        trim_main_title(title);
    }
    if title.ends_with(b", partial") {
        title.truncate(title.len() - 9);
        trim_main_title(title);
    }
    title.extend_from_slice(b" [");
    title.extend_from_slice(taxname);
    title.push(b']');
}

/// The title `CDeflineGenerator::GenerateDefline(fFastaFormat | fNoExpensiveOps)` makes of
/// the Bioseq's title descriptor `title` (`report::defline::generate_bioseq_defline`, with
/// the suffix adjustment of a protein with a source). `Err`: the Bioseqs whose title NCBI
/// makes from parts LOSAT does not model (an empty title made from the source; a source
/// without a name and a title ending in ` []`).
fn generated_title(source: &Source, set: &DeflineSet, title: &[u8]) -> Result<Vec<u8>, Stop> {
    let ids = &set.deflines[0].ids;
    let bioseq_ids = crate::report::db_subject::bioseq_ids(ids);
    if bioseq_title_from_source(title, &bioseq_ids) {
        return Err(Stop::Refused(format!(
            "-long_seqids of a sequence without a title (NCBI makes one from its source) is \
             not supported by LOSAT's {}",
            super::PROGRAM
        )));
    }
    let mut title = title.to_vec();
    if source.protein {
        // The first Source descriptor: the first definition line with a taxid above 0.
        let taxid = set
            .deflines
            .iter()
            .filter_map(|defline| defline.taxid)
            .find(|&taxid| taxid > 0);
        if let Some(taxid) = taxid {
            match source.tax_names(taxid)? {
                Some(names) => adjust_protein_title_suffix(&mut title, &names.scientific_name),
                None => {
                    let mut trimmed = title.clone();
                    trim_main_title(&mut trimmed);
                    if trimmed.ends_with(b" []") {
                        return Err(Stop::Refused(format!(
                            "-long_seqids of a protein whose title ends in \" []\" and whose \
                             taxid has no name is not supported by LOSAT's {}",
                            super::PROGRAM
                        )));
                    }
                }
            }
        }
    }
    let defline =
        generate_bioseq_defline(&title, source.protein, false, &bioseq_ids).map_err(|_| {
            Stop::Refused(format!(
                "a title whose encoding NCBI cannot tell is not supported by LOSAT's {}",
                super::PROGRAM
            ))
        })?;
    if defline.reads_past_end {
        return Err(Stop::Refused(format!(
            "a title on which NCBI's defline generator reads past its end is not supported \
             by LOSAT's {}",
            super::PROGRAM
        )));
    }
    Ok(defline.text)
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3295-3317
// ```c++
//             while ( count >= rem_line ) {
//                 x_WriteBuffer(ptr, rem_line);
//                 if ( !(current_state & eHardMask) ) {
//                     ptr += rem_line;
//                 }
//                 count -= rem_line;
//                 m_Out << '\n';
//                 rem_line = m_Width;
//             }
//             if ( count > 0 ) {
//                 x_WriteBuffer(ptr, count);
//                 rem_line -= count;
//             }
//             it.SetPos(new_pos);
//         }
//     }
//     if ( rem_line < m_Width ) {
//         m_Out << '\n';
//     }
// ```
/// The sequence lines of `CFastaOstream::x_WriteSequence`: lines of `width` letters, the
/// last one shorter, each ended by a line feed.
fn write_sequence(letters: &[u8], width: u32, out: &mut Vec<u8>) {
    for line in letters.chunks(width.max(1) as usize) {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
}

/// `ObjectId` of a `BL_ORD_ID` tag, for tests.
#[cfg(test)]
fn ordinal_id(n: i32) -> ParsedSeqId {
    ParsedSeqId::General(Dbtag {
        db: b"BL_ORD_ID".to_vec(),
        tag: ObjectId::Id(n),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqid::BlastDefline;

    fn defline(ids: Vec<ParsedSeqId>, title: Option<&str>) -> BlastDefline {
        BlastDefline {
            title: title.map(|t| t.as_bytes().to_vec()),
            ids,
            ..Default::default()
        }
    }

    fn text_ids(spec: &str) -> Vec<ParsedSeqId> {
        parse_ids(spec.as_bytes(), 0, &mut |_| {})
            .unwrap()
            .into_iter()
            .flatten()
            .collect()
    }

    #[test]
    fn the_title_joins_the_later_deflines_with_their_ids() {
        let set = DeflineSet {
            deflines: vec![
                defline(
                    text_ids("ref|NM_000001.1|"),
                    Some("alpha mRNA [Homo sapiens]"),
                ),
                defline(text_ids("gb|AB000001.1|"), Some("alpha cDNA")),
                defline(text_ids("emb|X00001.2|"), None),
                defline(vec![], Some("no ids")),
            ],
        };
        assert_eq!(
            bioseq_title(&set),
            b"alpha mRNA [Homo sapiens] >gb|AB000001.1| alpha cDNA >emb|X00001.2| no ids"
        );
        let configured = configure_defline_title(&bioseq_title(&set), false).unwrap();
        assert_eq!(
            configured,
            b"alpha mRNA [Homo sapiens] >AB000001.1 alpha cDNA >X00001.2 no ids"
        );
        let with_ctrl_a = configure_defline_title(&bioseq_title(&set), true).unwrap();
        assert_eq!(
            with_ctrl_a,
            b"alpha mRNA [Homo sapiens]\x01AB000001.1 alpha cDNA\x01X00001.2 no ids"
        );
        // A part that is not a Seq-id keeps its ` >`.
        assert_eq!(
            configure_defline_title(b"a >not an id >gi|5 x", true).unwrap(),
            b"a >not an id\x01gi|5 x"
        );
        assert_eq!(
            rejoin_defline_title(b"a >not an id >gi|5 x", true).unwrap(),
            b"a >not an id\x01gi|5 x"
        );
    }

    #[test]
    fn the_short_header_has_three_forms() {
        let header = |ids: Vec<ParsedSeqId>, range, minus| {
            let set = DeflineSet {
                deflines: vec![defline(ids, Some("title"))],
            };
            let config = Config {
                strand_minus: minus,
                ..Config::default()
            };
            let mut out = Vec::new();
            short_header(&set, range, &config, &mut out).unwrap();
            String::from_utf8(out).unwrap()
        };
        assert_eq!(header(vec![ordinal_id(3)], None, false), ">title\n");
        assert_eq!(
            header(vec![ordinal_id(3)], Some((2, 8)), true),
            ">:c9-3 title\n"
        );
        assert_eq!(
            header(text_ids("lcl|mynuc6"), None, false),
            ">mynuc6 title\n"
        );
        assert_eq!(
            header(text_ids("lcl|mynuc6"), Some((4, 19)), false),
            ">mynuc6:5-20 title\n"
        );
        assert_eq!(
            header(text_ids("gi|654321|gb|AB000002.3|"), None, false),
            ">AB000002.3 title\n"
        );
        assert_eq!(
            header(text_ids("gi|654321|gb|AB000002.3|"), Some((4, 19)), true),
            ">AB000002.3:c20-5 title\n"
        );
    }

    #[test]
    fn the_protein_suffix_is_rebuilt_from_the_organism() {
        let adjusted = |title: &str, taxname: &str| {
            let mut title = title.as_bytes().to_vec();
            adjust_protein_title_suffix(&mut title, taxname.as_bytes());
            String::from_utf8(title).unwrap()
        };
        assert_eq!(
            adjusted("alpha protein [Homo sapiens]", "Homo sapiens"),
            "alpha protein [Homo sapiens]"
        );
        assert_eq!(
            adjusted("alpha protein, partial [homo sapiens]", "Homo sapiens"),
            "alpha protein [Homo sapiens]"
        );
        assert_eq!(
            adjusted(
                "alpha protein (mitochondrion) [Homo sapiens]",
                "Homo sapiens"
            ),
            "alpha protein [Homo sapiens]"
        );
        assert_eq!(
            adjusted("alpha protein [Pan paniscus]", "Homo sapiens"),
            "alpha protein [Pan paniscus]"
        );
        assert_eq!(adjusted("alpha protein", "Homo sapiens"), "alpha protein");
    }

    #[test]
    fn sequence_lines_are_the_width_with_a_shorter_last() {
        let mut out = Vec::new();
        write_sequence(b"ACGTACGT", 3, &mut out);
        assert_eq!(out, b"ACG\nTAC\nGT\n");
        out.clear();
        write_sequence(b"ACGTAC", 3, &mut out);
        assert_eq!(out, b"ACG\nTAC\n");
    }

    #[test]
    fn split_keeps_empty_pieces() {
        assert_eq!(
            split_by_pattern(b"a >b > >c", b" >"),
            [&b"a"[..], b"b", b"", b"c"]
        );
        assert_eq!(split_by_pattern(b" >x", b" >"), [&b""[..], b"x"]);
        assert_eq!(split_by_pattern(b"", b" >"), [&b""[..]]);
    }
}
