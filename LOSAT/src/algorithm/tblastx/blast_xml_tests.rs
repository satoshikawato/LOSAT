use super::*;
use crate::stats::lookup_protein_params_ungapped;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:72-115: a query without a valid context has
// the -1 blocks and zero lengths.
#[test]
fn a_query_without_a_valid_context_has_the_unsearched_statistics() {
    let numbers = SubjectSetNumbers {
        total: 3_000,
        count: 3,
        link_total: 3_000,
        average: None,
        dbsize: None,
        searchsp: None,
        statistics: None,
    };
    assert_eq!(
        query_stats(&TblastxQueryStats::default(), &numbers),
        XmlStats::NONE
    );
}

// The length adjustment is the search's own function of the first valid context's length and
// the subject set (`-dbsize` replaces the length, not the number of sequences); the search
// space is the one the search kept.
#[test]
fn the_length_adjustment_is_computed_again_from_the_context_and_the_subject_set() {
    let karlin = lookup_protein_params_ungapped(crate::config::ScoringMatrix::Blosum62);
    for (dbsize, count) in [(None, 3), (Some(5_000_000), 3), (None, 1)] {
        let numbers = SubjectSetNumbers {
            total: 3_000,
            count,
            link_total: 3_000,
            average: None,
            dbsize,
            searchsp: None,
            statistics: None,
        };
        let lengths = try_compute_eff_lengths_tblastx(
            130,
            numbers.effective_length(),
            numbers.count,
            &karlin,
        )
        .expect("fits");
        let stats = TblastxQueryStats {
            karlin: Some(karlin),
            eff_searchsp: lengths.eff_searchsp,
            aa_len: 130,
            context0_search_space: 0,
            context0_aa_len: 0,
            seg_masks: Vec::new(),
        };
        let xml = query_stats(&stats, &numbers);
        assert_eq!(i64::from(xml.hsp_len), lengths.length_adjustment);
        assert!(xml.hsp_len > 0);
        assert_eq!(xml.eff_space, lengths.eff_searchsp);
        assert_eq!(
            (xml.kappa, xml.lambda, xml.entropy),
            (karlin.k, karlin.lambda, karlin.h)
        );
    }
}

// TBLASTX posts the Karlin-Altschul message only for a batch that was not searched; an empty
// record has its own message (blast_setup.hpp:190-198), and so has a query read past its earlier
// record (blast_setup_cxx.cpp:632-640), for which a batch not searched posts its message too
// (blast_aux_priv.cpp:94-101).
#[test]
fn the_messages_of_a_query_follow_its_batch_and_its_letters() {
    let query = FastaRecord::new("Query_1", b"q", b"ACGTACGTA");
    assert!(query_messages(&query, None, true, false).is_empty());
    assert_eq!(
        query_messages(&query, None, true, true),
        vec![EMPTY_DATABASE_MESSAGE.to_string()]
    );
    assert!(query_messages(&query, None, false, true)
        .iter()
        .all(|message| message != EMPTY_DATABASE_MESSAGE));
    assert_eq!(
        query_messages(&query, None, false, false),
        vec![INVALID_QUERY_MESSAGE.to_string()]
    );
    let empty = FastaRecord::new("Query_2", b"e", b"");
    assert_eq!(query_messages(&empty, None, true, false).len(), 1);
    assert_eq!(query_messages(&empty, None, false, false).len(), 2);
    let past_first = crate::blastinput::query_batch::PastFirstQuery {
        pos: 30,
        first_input: 0,
        first_length: 30,
    };
    assert_eq!(
        query_messages(&query, Some(past_first), true, false),
        vec!["CSeqVector_CI: cannot locate segment at 30".to_string()]
    );
    assert_eq!(
        query_messages(&query, Some(past_first), false, false),
        vec![
            "CSeqVector_CI: cannot locate segment at 30".to_string(),
            INVALID_QUERY_MESSAGE.to_string()
        ]
    );
    assert_eq!(
        setup_messages(&empty, None),
        vec![crate::algorithm::blastn::blast_engine::NO_DATA_MESSAGE.to_string()]
    );
}
