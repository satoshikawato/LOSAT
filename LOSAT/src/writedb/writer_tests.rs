use std::collections::BTreeSet;

use super::union_of_tax_ids;
use crate::seqid::BlastDefline;

/// A definition line with `taxid` and `links`.
fn line(taxid: Option<i32>, links: Option<&[i32]>) -> BlastDefline {
    BlastDefline {
        taxid,
        links: links.map(<[i32]>::to_vec),
        ..BlastDefline::default()
    }
}

/// `m_TaxIds` of a set gathered in one vector is the union of each line's `GetTaxIds`
/// (`BlastDefline::tax_ids`): a taxid alone, a taxid 0 alone or beside links, links that hold 0
/// (beside other values, or only 0s), repeated values, and two or three lines whose 0 stays only
/// where a line has nothing else.
#[test]
fn the_tax_ids_of_a_set_are_the_union_of_each_line_s() {
    let cases: Vec<Vec<BlastDefline>> = vec![
        vec![],
        vec![line(None, None)],
        vec![line(Some(9606), None)],
        vec![line(Some(0), None)],
        vec![line(Some(0), Some(&[562, 7227]))],
        vec![line(Some(9606), Some(&[0, 562]))],
        vec![line(None, Some(&[0, 0]))],
        vec![line(Some(0), Some(&[0]))],
        vec![line(Some(562), Some(&[562, 562]))],
        vec![line(None, Some(&[]))],
        vec![line(Some(0), None), line(Some(9606), None)],
        vec![line(Some(0), Some(&[10090])), line(Some(0), None)],
        vec![
            line(Some(7227), Some(&[0, 562])),
            line(Some(562), None),
            line(None, Some(&[0, 9606])),
        ],
    ];
    for deflines in cases {
        let expected: BTreeSet<i32> = deflines.iter().flat_map(BlastDefline::tax_ids).collect();
        assert_eq!(
            union_of_tax_ids(&deflines),
            expected.into_iter().collect::<Vec<_>>(),
            "{deflines:?}"
        );
    }
}
