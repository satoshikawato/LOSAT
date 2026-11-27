//! Unit tests of the TSE feature index of `asn_defline`: the Seq-id interner's chains of entries
//! with one hash, and the picks NCBI's scope may take among candidate Bioseqs.

use super::*;

/// Two Local Seq-ids whose FxHashes are equal (rustc-hash 2.1.3): the interner keeps each once and
/// finds both again through the chain of entries with that hash.
#[test]
fn interner_keeps_seq_ids_with_one_hash_apart() {
    let a = ParsedSeqId::Local(ObjectId::Str(b"probe-local-id-1".to_vec()));
    let b = ParsedSeqId::Local(ObjectId::Str(vec![
        0xf4, 0x1a, 0xbf, 0xab, 0xcf, 0x84, 0x0b, 0x06, 0xe7, 0x09, 0xbc, 0xe4, 0xc3, 0xcd, 0x4a,
        0x58,
    ]));
    let hash = |id: &ParsedSeqId| {
        let mut hasher = FxHasher::default();
        id.hash(&mut hasher);
        hasher.finish()
    };
    assert_eq!(
        hash(&a),
        hash(&b),
        "the pair no longer collides: pick another"
    );
    assert_ne!(a, b);
    let mut interner = SeqIdInterner::default();
    let ra = interner.intern(&a);
    let rb = interner.intern(&b);
    assert_ne!(ra, rb);
    assert_eq!(
        interner.same_hash[rb as usize], ra,
        "b's chain goes on to a"
    );
    assert_eq!(interner.same_hash[ra as usize], NO_ID);
    // Found again past the last entry, at the head of the chain and further down it.
    for (id, r) in [(&a, ra), (&b, rb), (&a, ra), (&a, ra), (&b, rb)] {
        assert_eq!(interner.intern(id), r);
    }
    assert_eq!(interner.ids.len(), 2, "each Seq-id is held once");
    assert_eq!(interner.id(ra), &a);
    assert_eq!(interner.id(rb), &b);
}

fn cand(bio: u32, tree: u8, version: Option<i32>) -> Candidate {
    Candidate { bio, tree, version }
}

fn picks(found: Found) -> Vec<Option<u32>> {
    found_picks(found).into_vec()
}

/// PDB Seq-ids have no better version (`IsBetterVersion` is false): two Bioseqs of one tree whose
/// PDB Seq-ids the scope finds alike are both picks; text Seq-ids of one tree keep the best
/// version; candidates of two trees are all picks.
#[test]
fn picks_follow_the_scope_match_loop() {
    assert_eq!(picks(pick_candidate(&[])), vec![None]);
    assert_eq!(picks(pick_candidate(&[cand(3, 7, None)])), vec![Some(3)]);
    assert_eq!(
        picks(pick_candidate(&[cand(1, 7, None), cand(2, 7, None)])),
        vec![Some(1), Some(2)]
    );
    assert_eq!(
        picks(pick_candidate(&[cand(1, 2, Some(1)), cand(2, 2, Some(2))])),
        vec![Some(2)]
    );
    assert_eq!(
        picks(pick_candidate(&[cand(1, 2, Some(1)), cand(2, 3, Some(2))])),
        vec![Some(1), Some(2)]
    );
    // The same Bioseq twice is one pick.
    assert_eq!(
        picks(pick_candidate(&[cand(4, 2, Some(1)), cand(4, 3, None)])),
        vec![Some(4)]
    );
}

/// Name matches that NCBI may or may not take (`accession_name_match`) add the picks of every
/// subset of their trees: none (here no Bioseq) and each.
#[test]
fn unsure_name_matches_add_both_outcomes() {
    assert_eq!(
        picks(pick_candidate_unsure(&[], &[cand(5, 9, None)])),
        vec![None, Some(5)]
    );
    assert_eq!(
        picks(pick_candidate_unsure(
            &[cand(1, 2, Some(1))],
            &[cand(5, 9, None)]
        )),
        vec![Some(1), Some(5)]
    );
    // Without unsure matches it is `pick_candidate`.
    assert_eq!(
        picks(pick_candidate_unsure(&[cand(1, 2, None)], &[])),
        vec![Some(1)]
    );
}
