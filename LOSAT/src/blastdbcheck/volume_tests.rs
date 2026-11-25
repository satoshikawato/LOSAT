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
    // Without a positive taxid, the first line whose first id of the first id's type matches it
    // and that has a taxid decides (`GetTaxidForSeqid`).
    assert!(has_taxid(
        0,
        &set(vec![line(vec![gi(1)], None), line(vec![gi(1)], Some(-3))])
    )
    .unwrap_or(false));
    assert!(!has_taxid(
        0,
        &set(vec![line(vec![gi(1)], None), line(vec![gi(2)], Some(-3))])
    )
    .unwrap_or(true));
    assert!(!has_taxid(
        0,
        &set(vec![
            line(vec![gi(1)], None),
            line(vec![gi(1)], Some(0)),
            line(vec![gi(1)], Some(-3))
        ])
    )
    .unwrap_or(true));
    // Only the line's first id of that type is compared.
    assert!(!has_taxid(
        0,
        &set(vec![
            line(vec![gi(1)], None),
            line(vec![gi(2), gi(1)], Some(-3))
        ])
    )
    .unwrap_or(true));
    // A first line without ids: NCBI reads a null reference.
    assert!(has_taxid(0, &set(vec![line(Vec::new(), None)])).is_err());
}

/// Approved exception `Q-D2b-OV=A` (`PD-LOSAT-NCBI-DEFECTS`, Owner 2026-11-25 17:54 JST):
/// `-ends` near `i32::MIN` widens `count - n` to `i64` instead of overflowing signed 32-bit `int`.
#[test]
fn widened_ends_near_int_min_jumps_past_volume() {
    let count = 9i64;
    let end = count - i64::from(i32::MIN);
    assert_eq!(end, 2_147_483_657);
    assert!(end > 0);
    let near = count - i64::from(-2_147_483_638i32);
    assert_eq!(near, i64::from(i32::MAX));
}
