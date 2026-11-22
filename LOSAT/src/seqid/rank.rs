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
/// `CSeq_id::Score` (= `WorstRank`, `TextScore`): lower is better.
pub fn text_score(id: &ParsedSeqId) -> i32 {
    adjust_score(id, base_text_score(id))
}

/// `CSeq_id::BestRank` (`BestRankScore`): lower is better.
pub fn best_rank_score(id: &ParsedSeqId) -> i32 {
    adjust_score(id, base_best_rank_score(id))
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3126-3153
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
/// `FindBestChoice` of `ids` with `score`: the first id with the lowest score; `None` for no ids.
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
