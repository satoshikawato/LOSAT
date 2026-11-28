//! Unit tests of `seqid::relation`: `CSeq_id::Match` of the Seq-ids that `x_FillLocList` compares,
//! a mask's and an alignment row's (line RG-KB item 12d).

use super::*;
use crate::seqdb::Oid;
use crate::seqid::accguide::parse_flags;
use crate::seqid::parse::seq_id_from_text;

/// The Seq-id of the FASTA id `text`.
fn id(text: &str) -> ParsedSeqId {
    seq_id_from_text(text.as_bytes(), parse_flags::ANY_RAW, &mut |_| {}).expect("a Seq-id")
}

// The query `gb|KCW91010|` of a TBLASTX search with `-parse_deflines` and the database sequence
// `gb|KCW91010.1|` (line RG-KB item 12c, case `ctl_rev_2`): NCBI 2.17.0 placed the query's masks
// through the subject's row, so the two ids match although NCBI's scope does not find the query
// for the subject.
#[test]
fn a_text_id_matches_without_a_version_and_across_the_genbank_kinds() {
    let matches = |a: &str, b: &str| super::matches(&id(a), &id(b));
    assert!(matches("gb|KCW91010|", "gb|KCW91010.1|"));
    assert!(matches("gb|KCW91010.1|", "gb|KCW91010|"));
    assert!(matches("gb|kcw91010.1|", "gb|KCW91010.1|"));
    assert!(matches("emb|KCW91010|", "gb|KCW91010.1|"));
    assert!(!matches("gb|KCW91010.1|", "gb|KCW91010.2|"));
    assert!(!matches("gb|KCW91011|", "gb|KCW91010.1|"));
    assert!(!matches("ref|NC_000913.3|", "gb|NC_000913.3|"));
    // NCBI's scope does not find the query `gb|KCW91010|` for the subject's id: the masks are
    // placed by `matches`, not by the scope's finding.
    assert!(!scope_finds(&id("gb|KCW91010|"), &id("gb|KCW91010.1|")));
}

// The generated ids of a run without `-parse_deflines` (`Query_N`, `Subject_N`) are local ids;
// a database sequence without parsed Seq-ids is `gnl|BL_ORD_ID|n`.
#[test]
fn reported_ids_match_as_their_seq_ids() {
    let local = |text| SeqId::Local(text);
    let oid = |n| SeqId::BlOrdId(Oid::of_unrestricted(n));
    let lcl = id("lcl|Query_1");
    assert!(row_matches(local("Query_1"), local("Query_1")));
    assert!(!row_matches(local("Query_1"), local("Subject_1")));
    assert!(row_matches(local("Query_1"), SeqId::Db(&lcl)));
    assert!(row_matches(oid(5), oid(5)));
    assert!(!row_matches(oid(5), oid(6)));
    assert!(!row_matches(local("Query_1"), oid(0)));
    let general = id("gnl|BL_ORD_ID|5");
    assert!(row_matches(SeqId::Db(&general), oid(5)));
    let text = id("gb|KCW91010.1|");
    assert!(!row_matches(SeqId::Db(&text), oid(5)));
    assert!(!row_matches(SeqId::Db(&text), local("Query_1")));
}

#[test]
fn row_pairs_hold_the_pairs_whose_ids_match() {
    let queries = [id("gb|KCW91010|"), id("gb|KCW91011|")];
    let subjects = [id("gb|KCW91010.1|"), id("gb|KCW91012.1|")];
    let pairs = RowPairs::of(
        [(0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (0, 0)],
        |q| SeqId::Db(&queries[q as usize]),
        |s| SeqId::Db(&subjects[s as usize]),
    );
    assert!(pairs.matches(0, 0));
    assert!(!pairs.matches(0, 1));
    assert!(!pairs.matches(1, 0));
    assert!(!pairs.matches(1, 1));
    assert!(!RowPairs::default().matches(0, 0));
}
