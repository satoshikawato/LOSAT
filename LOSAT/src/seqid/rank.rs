//! The ranks by which NCBI chooses one Seq-id of a list (`FindBestChoice`): `Score` (also
//! `WorstRank`, the `TextScore`), which the formatters use for the id they show, and
//! `BestRank`, which picks the Seq-align id of a database sequence (`defline_set`).
//!
//! The rank of a `gi` id depends on two environment variables of NCBI's
//! (`SEQ_ID_PREFER_ACCESSION_OVER_GI`, `SEQ_ID_AVOID_GI`, Seq_id.cpp:3393-3413), which make it
//! the worst; these functions give the ranks with both unset, NCBI's default. LOSAT does not
//! refuse the two variables yet: no output reads these ranks before line DBID ID2a, which
//! adds them to `blastinput::ncbi_environment`'s refused variables with the first one (the
//! registry spellings are refused already).

use super::model::ParsedSeqId;

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3039-3077
// ```c++
// int CSeq_id::BaseTextScore(void) const
// {
//     switch (Which()) {
//     // Accession and accession-like ids - only one can be present in a bioseq's list of ids,
//     // the order is not important.
//     case e_Other:
//     case e_Swissprot:
//     case e_Pir:
//     case e_Pdb:
//     case e_Genbank:
//     case e_Embl:
//     case e_Ddbj:
//     case e_Tpg:
//     case e_Tpe:
//     case e_Tpd:                 return 10;
//
//     // Second group of mutually exclusive ids, any order can be used.
//     case e_Gpipe:
//     case e_Named_annot_track:
//     case e_Prf:                 return 20;
//     case e_Patent:              return 50;
//
//     // "local" < "general" < "gi"
//     case e_Local:               return 50;
//     case e_General:             return 60;
//     case e_Gi:                  return PreferAccessionOverGi() ? kMaxScore + 1 : 70;
//
//     // All other ids rank just above "not-set" except the obsolete ones listed below.
//     case e_not_set:             return 100;
//
//     // Obsolete ids, lowest rank, any order is OK.
//     case e_Giim:
//     case e_Gibbmt:
//     case e_Gibbsq:              return 1000;
//
//     // All other ids should go just above "not-set".
//     default:                    return 90;
//     }
// }
// ```
fn base_text_score(id: &ParsedSeqId) -> i32 {
    match id {
        ParsedSeqId::Other(_)
        | ParsedSeqId::Swissprot(_)
        | ParsedSeqId::Pir(_)
        | ParsedSeqId::Pdb(_)
        | ParsedSeqId::Genbank(_)
        | ParsedSeqId::Embl(_)
        | ParsedSeqId::Ddbj(_)
        | ParsedSeqId::Tpg(_)
        | ParsedSeqId::Tpe(_)
        | ParsedSeqId::Tpd(_) => 10,
        ParsedSeqId::Gpipe(_) | ParsedSeqId::NamedAnnotTrack(_) | ParsedSeqId::Prf(_) => 20,
        ParsedSeqId::Patent(_) => 50,
        ParsedSeqId::Local(_) => 50,
        ParsedSeqId::General(_) => 60,
        ParsedSeqId::Gi(_) => 70,
        ParsedSeqId::Giim(_) | ParsedSeqId::Gibbmt(_) | ParsedSeqId::Gibbsq(_) => 1000,
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3080-3093
// ```c++
// int CSeq_id::BaseBestRankScore(void) const
// {
//     switch (Which()) {
//     case e_not_set:                               return 83;
//     case e_General: case e_Local:                 return 80;
//     case e_Gibbsq: case e_Gibbmt: case e_Giim:    return 70;
//     case e_Named_annot_track:                     return 69;
//     case e_Gpipe:                                 return 68;
//     case e_Patent:                                return 67;
//     case e_Other:                                 return 65;
//     case e_Gi: return PreferAccessionOverGi() ? kMaxScore + 1 : 51;
//     default:                                      return 60;
//     }
// }
// ```
fn base_best_rank_score(id: &ParsedSeqId) -> i32 {
    match id {
        ParsedSeqId::General(_) | ParsedSeqId::Local(_) => 80,
        ParsedSeqId::Gibbsq(_) | ParsedSeqId::Gibbmt(_) | ParsedSeqId::Giim(_) => 70,
        ParsedSeqId::NamedAnnotTrack(_) => 69,
        ParsedSeqId::Gpipe(_) => 68,
        ParsedSeqId::Patent(_) => 67,
        ParsedSeqId::Other(_) => 65,
        ParsedSeqId::Gi(_) => 51,
        _ => 60,
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3195-3221
// ```c++
// int CSeq_id::AdjustScore(int base_score, TAdjustScoreFlags flags) const
// {
//     int score = base_score * 10;
//     if ( IsGeneral() ) {
//         const string& db = GetGeneral().GetDb();
//         if ( db == "TRACE" ) {
//             // prefer "ti" over "TRACE"
//             score += 5;
//         }
//     }
//     else if ( const CTextseq_id* text_id = GetTextseq_Id() ) {
//         if ( !text_id->IsSetVersion() ) {
//             score += 4;
//         }
//         if ( !text_id->IsSetAccession() ) {
//             if ((flags & fRequireAccessions) == 0) {
//                 score += 3; // still penalize somewhat
//             } else {
//                 score = kMax_Int;
//             }
//         }
//         if ( !text_id->IsSetName() ) {
//             score += 2;
//         }
//     }
//     return score;
// }
// ```
/// `AdjustScore` without `fRequireAccessions` (the ranks below do not set it).
fn adjust_score(id: &ParsedSeqId, base_score: i32) -> i32 {
    let mut score = base_score * 10;
    if let ParsedSeqId::General(dbtag) = id {
        if dbtag.db == b"TRACE" {
            score += 5;
        }
    } else if let Some(text_id) = id.textseq() {
        if text_id.version.is_none() {
            score += 4;
        }
        if text_id.accession.is_none() {
            score += 3;
        }
        if text_id.name.is_none() {
            score += 2;
        }
    }
    score
}

// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:753-779
// ```c++
//     int TextScore     (void) const { return AdjustScore(BaseTextScore()); }
//     int BestRankScore (void) const { return AdjustScore(BaseBestRankScore()); }
//     int WorstRankScore(void) const { return TextScore(); }
//     ...
//     static int Score(const CRef<CSeq_id>& id)
//         { return id ? id->TextScore() : kMax_Int; }
//     ...
//     static int BestRank(const CRef<CSeq_id>& id)
//         { return id ? id->BestRankScore() : kMax_Int; }
//     static int WorstRank(const CRef<CSeq_id>& id)
//         { return Score(id); }
// ```
/// `CSeq_id::Score` (= `WorstRank`, `TextScore`): lower is better.
pub fn text_score(id: &ParsedSeqId) -> i32 {
    adjust_score(id, base_text_score(id))
}

/// `CSeq_id::BestRank` (`BestRankScore`): lower is better.
pub fn best_rank_score(id: &ParsedSeqId) -> i32 {
    adjust_score(id, base_best_rank_score(id))
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3126-3153
// ```c++
// int CSeq_id::BaseFastaAAScore(void) const
// {
//     switch (Which()) {
//     case e_not_set: case e_Giim:   return 255;
//     case e_Local:                  return 230;
//     case e_Gi: return PreferAccessionOverGi() ? kMaxScore + 1 : 120;
//     case e_General:
//         {
//         const string& db = GetGeneral().GetDb();
//         if (db.compare("TMSMART") == 0 ||
//             db.compare("BankIt") == 0 ||
//             db.compare("NCBIFILE") == 0)
//             return 240;
//         else
//             return 90;
//         }
//     case e_Patent:                 return 80;
//     case e_Prf:                    return 70;
//     case e_Pdb:                    return 50;
// //  see SQD-4175 ticket for priorities
//     case e_Gibbsq:                 return 42;
//     case e_Gibbmt:                 return 41;
//     case e_Genbank:                return 40;
//     case e_Pir:                    return 30;
//     case e_Swissprot:              return 20;
//     case e_Other:                  return 15;
//     default:                       return 60; // [third party] GB/EMBL/DDBJ
//     }
// }
// ```
fn base_fasta_aa_score(id: &ParsedSeqId) -> i32 {
    match id {
        ParsedSeqId::Giim(_) => 255,
        ParsedSeqId::Local(_) => 230,
        ParsedSeqId::Gi(_) => 120,
        ParsedSeqId::General(dbtag) => {
            if [&b"TMSMART"[..], b"BankIt", b"NCBIFILE"].contains(&dbtag.db.as_slice()) {
                240
            } else {
                90
            }
        }
        ParsedSeqId::Patent(_) => 80,
        ParsedSeqId::Prf(_) => 70,
        ParsedSeqId::Pdb(_) => 50,
        ParsedSeqId::Gibbsq(_) => 42,
        ParsedSeqId::Gibbmt(_) => 41,
        ParsedSeqId::Genbank(_) => 40,
        ParsedSeqId::Pir(_) => 30,
        ParsedSeqId::Swissprot(_) => 20,
        ParsedSeqId::Other(_) => 15,
        _ => 60,
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3096-3124
// ```c++
// int CSeq_id::BaseFastaNAScore(void) const
// {
//     switch (Which()) {
//         // these few are bogus, at least for nucleotide sequences
//     case e_not_set: case e_Giim:
//     case e_Pir: case e_Swissprot: case e_Prf:  return 255;
//     case e_Local:                              return 230;
//     case e_Gi: return PreferAccessionOverGi() ? kMaxScore + 1 : 120;
//     case e_General:
//         {
//         const string& db = GetGeneral().GetDb();
//         if (db.compare("TMSMART") == 0 ||
//             db.compare("BankIt") == 0 ||
//             db.compare("NCBIFILE") == 0 )
//             return 240;
//         else
//             return 100;
//         }
//     case e_Patent:                 return 90;
//     case e_Pdb:                    return 80;
// //  see SQD-4175 ticket for priorities
//     case e_Gibbsq:                 return 72;
//     case e_Gibbmt:                 return 71;
//     case e_Genbank:                return 70;
//     case e_Other:                  return 15;
//     default: /* [third party] GB/EMBL/DDBJ */  return 20;
//     }
// }
// ```
fn base_fasta_na_score(id: &ParsedSeqId) -> i32 {
    match id {
        ParsedSeqId::Giim(_)
        | ParsedSeqId::Pir(_)
        | ParsedSeqId::Swissprot(_)
        | ParsedSeqId::Prf(_) => 255,
        ParsedSeqId::Local(_) => 230,
        ParsedSeqId::Gi(_) => 120,
        ParsedSeqId::General(dbtag) => {
            if [&b"TMSMART"[..], b"BankIt", b"NCBIFILE"].contains(&dbtag.db.as_slice()) {
                240
            } else {
                100
            }
        }
        ParsedSeqId::Patent(_) => 90,
        ParsedSeqId::Pdb(_) => 80,
        ParsedSeqId::Gibbsq(_) => 72,
        ParsedSeqId::Gibbmt(_) => 71,
        ParsedSeqId::Genbank(_) => 70,
        ParsedSeqId::Other(_) => 15,
        _ => 20,
    }
}

/// `CSeq_id::FastaAARank` (`FastaAAScore`): the id `CFastaOstream` writes first for a
/// protein; lower is better.
pub fn fasta_aa_score(id: &ParsedSeqId) -> i32 {
    adjust_score(id, base_fasta_aa_score(id))
}

/// `CSeq_id::FastaNARank` (`FastaNAScore`): the same for a nucleotide.
pub fn fasta_na_score(id: &ParsedSeqId) -> i32 {
    adjust_score(id, base_fasta_na_score(id))
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbiutil.hpp:224-258
// ```c++
//     void operator() (const T& x)
//     {
//         int score = m_Func(x);
//         if (score < m_Score) {
//             m_Value = x;
//             m_Score = score;
//         }
//     }
//     ...
// template <typename C, typename F>
// inline
// typename C::value_type
// FindBestChoice(const C& container, F score_func)
// {
//     typedef typename C::value_type T;
//     CBestChoiceTracker<T, F> tracker(score_func);
//     ITERATE (typename C, it, container) {
//         tracker(*it);
//     }
//     return tracker.GetBestChoice();
// }
// ```
/// `FindBestChoice(ids, score)`: the first id with the lowest score; `None` for no ids.
pub fn find_best<'a>(
    ids: impl IntoIterator<Item = &'a ParsedSeqId>,
    score: impl Fn(&ParsedSeqId) -> i32,
) -> Option<&'a ParsedSeqId> {
    let mut best = None;
    let mut best_score = i32::MAX;
    for id in ids {
        let value = score(id);
        if value < best_score {
            best = Some(id);
            best_score = value;
        }
    }
    best
}
