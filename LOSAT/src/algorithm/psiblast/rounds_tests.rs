use super::*;

fn set(keys: &[u32]) -> InclusionSet {
    keys.iter().map(|&key| SetKey::Index(key)).collect()
}

// psiblast_iteration.cpp:68-96 (cited at `IterationState`): no convergence before two
// rounds, an empty set after a non-empty one, and every current id among the previous ones
// (two empty sets too).
#[test]
fn has_converged_follows_ncbis_three_rules() {
    let mut state = IterationState::new(0);
    assert!(!state.has_converged());
    state.advance(set(&[]));
    assert!(!state.has_converged(), "one round: never converged");
    state.advance(set(&[]));
    assert!(state.has_converged(), "two empty sets converge");
    let mut state = IterationState::new(0);
    state.advance(set(&[1, 2]));
    assert!(!state.has_converged());
    state.advance(set(&[]));
    assert!(state.has_converged(), "an empty set after a non-empty one");
    let mut state = IterationState::new(0);
    state.advance(set(&[1, 2]));
    state.advance(set(&[2]));
    assert!(state.has_converged(), "a subset of the previous set");
    let mut state = IterationState::new(0);
    state.advance(set(&[1, 2]));
    state.advance(set(&[2, 3]));
    assert!(!state.has_converged(), "a new id");
    assert!(state.go(), "no limit");
}

// psiblast_iteration.cpp:58-65, 113-129: the number of rounds, the round's number and the
// previous set that the report reads.
#[test]
fn the_rounds_run_out_and_number_from_one() {
    let mut state = IterationState::new(2);
    assert!(state.go());
    assert_eq!(state.number(), 1);
    assert!(state.previously_found().is_empty());
    state.advance(set(&[7]));
    assert_eq!(state.number(), 2);
    assert_eq!(state.previously_found(), &set(&[7]));
    assert!(state.go());
    state.advance(set(&[7, 8]));
    assert!(!state.has_more_iterations());
    assert!(!state.go());
    assert!(!state.has_converged());
    let mut state = IterationState::new(2);
    state.advance(set(&[7]));
    state.advance(set(&[7]));
    assert!(state.has_converged(), "converged on the last allowed round");
}

// psiblast_app.cpp:296-310 (cited at `blastp::BlastpReport::resolve_subjects`): a round makes a
// PSSM, which reads its subjects' letters from NCBI's scope, only where another round can follow
// it or the last PSSM is saved (review 3 of line RG-KB item 12, F1).
#[test]
fn a_round_may_make_a_pssm_only_where_a_round_follows_or_the_last_is_saved() {
    let io = |total: u32, done: u32, save_last: bool| {
        let options = PsiblastPssmOptions {
            rounds: Some(total),
            save_pssm_after_last_round: save_last,
            ..PsiblastPssmOptions::default()
        };
        let plan = PsiRoundPlan::of(&options, None).expect("rounds");
        let mut state = IterationState::new(total);
        for round in 0..done {
            state.advance(set(&[round]));
        }
        PsiRoundIo {
            state,
            save_last: plan.save_last,
            request: plan.request,
            query_id: None,
            outcome: std::sync::Mutex::new(None),
        }
    };
    assert!(!io(1, 0, false).may_make_next(), "one round");
    assert!(
        io(1, 0, true).may_make_next(),
        "one round whose PSSM is saved"
    );
    assert!(io(2, 0, false).may_make_next(), "round 1 of 2");
    assert!(!io(2, 1, false).may_make_next(), "round 2 of 2");
    assert!(io(0, 3, false).may_make_next(), "until convergence");
}

fn id(text: &str) -> SetKey {
    use crate::seqid::accguide::parse_flags::{ANY_LOCAL, PARTIAL_OK, RAW_TEXT};
    use crate::seqid::parse::parse_ids;
    let ids = parse_ids(
        text.as_bytes(),
        PARTIAL_OK | ANY_LOCAL | RAW_TEXT,
        &mut |_| {},
    )
    .unwrap();
    assert_eq!(ids.len(), 1, "{text}");
    SetKey::Id(ids.into_iter().next().unwrap().unwrap())
}

#[test]
fn keys_of_records_and_ids_differ() {
    let ids: InclusionSet = [SetKey::Index(3), id("gb|X1.1|")].into_iter().collect();
    assert!(ids.contains(&SetKey::Index(3)));
    assert!(ids.contains(&id("gb|X1.1|")));
    assert!(!ids.contains(&id("gb|X1|")));
    assert!(!ids.contains(&id("lcl|3")));
}

// Line DBID pssmrows (weakbatch DECISIONS 43 (4)): the keys of a database's subjects are their
// Seq-align ids as `ParsedSubject` reads them; psiblast's PSSM puts following hits whose ids
// `Match` in one row and counts its rows by their handles (psi_pssm_input.cpp
// x_CountAndSelectQualifyingAlignments and x_ExtractAlignmentData, cited at
// `pssm::input_msa`).
#[test]
fn a_pssm_row_joins_following_hits_whose_seq_align_ids_match() {
    use crate::seqdb::{SeqDb, SeqType};
    // tests/db_search/ids/pweak_v5: 7 `gb|FFF00043|`, 8 `gb|FFF00043.1|`, 4 `gb|FFF00041.1|`.
    let base =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/pweak_v5");
    let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
    let key = |ordinal: u32| {
        let subject =
            crate::report::db_subject::ParsedSubject::new(db.defline_set(ordinal).unwrap())
                .unwrap();
        SetKey::Id(subject.alignment_id().clone())
    };
    // Two handles whose ids `Match`: two keys, one row where they follow each other.
    assert_ne!(key(8), key(7));
    assert!(key(8).matches(&key(7)) && key(7).matches(&key(8)));
    assert!(!key(8).matches(&key(4)));
    SetKey::check_handles(&[key(8), key(7), key(4)]).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:310-320 (cited at
// `pssm::input_msa`): following Seq-aligns whose ids `Match` share a row.
#[test]
fn keys_match_as_their_seq_align_ids() {
    assert!(id("gb|X1|").matches(&id("gb|X1.1|")));
    assert!(id("gb|X1.1|").matches(&id("emb|x1.1|")));
    assert!(!id("gb|X1.1|").matches(&id("gb|X1.2|")));
    assert!(id("sp||NAMEY").matches(&id("sp|P00042.1|NAMEY")));
    assert!(SetKey::Index(5).matches(&SetKey::Index(5)));
    assert!(!SetKey::Index(5).matches(&SetKey::Index(6)));
    assert!(SetKey::Index(5).matches(&id("gnl|BL_ORD_ID|5")));
    assert!(!id("gnl|BL_ORD_ID|5").matches(&SetKey::Index(6)));
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1064-1095 (cited at
// `seqid::defline_set::handle_superset_key`): a handle compares a text id's accession
// without case.
#[test]
fn keys_that_ncbis_handles_can_take_for_one_are_refused() {
    let keys = [
        id("gb|X1|"),
        id("gb|X1.1|"),
        id("gb|X1.1|"),
        id("sp||NAMEY"),
    ];
    SetKey::check_handles(&keys).unwrap();
    SetKey::check_handles(&[id("pdb|1ABC|A"), id("pdb|1ABC|a"), id("pdb|1ABC|B")]).unwrap();
    for (a, b) in [
        ("gb|ABC12.1|", "gb|abc12.1|"),
        ("gb|ABC12.1|", "emb|ABC12.1|"),
        ("lcl|seqA", "lcl|SEQA"),
        ("gnl|MyDb|tagX", "gnl|MYDB|TAGX"),
    ] {
        let error = SetKey::check_handles(&[id(a), id("gb|X1.1|"), id(b)]).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!("the Seq-align ids {a} and {b} of two database sequences that psiblast compares (in its inclusion sets, the rows of its PSSM or its two-part description table), which NCBI BLAST+'s Seq-id handles can take for one sequence (ids that differ only in letter case, one accession under gb, emb and dbj, and other ids that LOSAT does not tell apart as NCBI's handles do), are not supported by LOSAT's PSIBLAST")
        );
    }
    let error = SetKey::check_handles(&[SetKey::Index(5), id("gnl|BL_ORD_ID|5")]).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("gnl|BL_ORD_ID|5 (a sequence without parsed Seq-ids) and gnl|BL_ORD_ID|5 of"),
        "{error}"
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2926-2938 (cited at
// `HandleIndex`): a giim id is never one handle with an equal one.
#[test]
fn a_giim_key_is_refused() {
    let error = SetKey::check_handles(&[id("gb|X1.1|"), id("gim|123")]).unwrap_err();
    assert!(error.to_string().contains("gim|123"), "{error}");
    assert!(error.to_string().contains("a giim id"), "{error}");
}

// Review 1 of line DBID pssmrows, M1: a key outside the index is unsure only where it shares
// its superset key with an indexed key that differs.
#[test]
fn a_key_is_unsure_against_a_previous_set_only_where_the_handles_can_be_one() {
    let previous = [id("gb|ABC12.1|"), id("gb|X1.1|")];
    let index = HandleIndex::new(&previous).unwrap();
    assert_eq!(index.unsure(&id("gb|abc12.1|")), Some(&previous[0]));
    assert_eq!(index.unsure(&id("emb|ABC12.1|")), Some(&previous[0]));
    assert_eq!(index.unsure(&id("gb|ABC12.1|")), None);
    assert_eq!(index.unsure(&id("gb|X1|")), None);
    assert_eq!(index.unsure(&id("gb|Y1.1|")), None);
}
