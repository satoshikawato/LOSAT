use super::*;
use crate::seqid::defline_set::BlastDefline;
use crate::seqid::ParsedSeqId;

fn line(ids: Vec<ParsedSeqId>, taxid: Option<i32>) -> BlastDefline {
    BlastDefline {
        ids,
        taxid,
        ..BlastDefline::default()
    }
}

fn set(lines: Vec<BlastDefline>) -> DeflineSet {
    DeflineSet::from(lines)
}

/// A positive taxid of any line makes the Bioseq's source; else the first line's taxid decides
/// (probes ca4/ca5: the second line's 9606 passes; taxids_p5 OID 6: none fails).
#[test]
fn a_taxid_is_any_positive_one_else_the_first_lines() {
    let gi = |n| ParsedSeqId::Gi(n);
    assert!(has_taxid(
        0,
        &set(vec![line(vec![gi(1)], None), line(vec![gi(2)], Some(9606))])
    )
    .unwrap_or(false));
    assert!(!has_taxid(0, &set(vec![line(vec![gi(1)], None)])).unwrap_or(true));
    assert!(!has_taxid(0, &set(vec![line(vec![gi(1)], Some(0))])).unwrap_or(true));
    assert!(has_taxid(0, &set(vec![line(vec![gi(1)], Some(-3))])).unwrap_or(false));
    // Another line of a negative taxid with an id of the first id's type: CSeq_id::Match decides.
    assert!(has_taxid(
        0,
        &set(vec![line(vec![gi(1)], None), line(vec![gi(1)], Some(-3))])
    )
    .is_err());
    // A first line without ids: NCBI reads a null reference.
    assert!(has_taxid(0, &set(vec![line(Vec::new(), None)])).is_err());
}
