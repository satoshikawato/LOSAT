use super::*;

/// A record written by `encode` reads back with the same fields.
#[test]
fn capture_record_round_trips() {
    let q = [3u8, 5, 7, 1];
    let s = [2u8, 4, 6, 8, 10];
    let m: Vec<i32> = (0..16 * 16).map(|i| (i % 7) - 3).collect();
    let rows = RowSeq::Bytes {
        data: &q,
        base: 4,
        step: -1,
    };
    let cols = ColSeq {
        data: &s,
        base: 4,
        step: -1,
        zero_from: 4,
    };
    let scores = Scores::Flat { data: &m, n: 16 };
    let res = XdropResult {
        a_offset: 3,
        b_offset: 2,
        score: 17,
        fence_hit: false,
        cells: 99,
    };
    let ops = [(3u8, 2u32), (0, 1)];
    let buf = encode(
        &rows,
        &cols,
        &scores,
        4,
        4,
        11,
        1,
        38,
        true,
        true,
        Some(&res),
        &ops,
    )
    .unwrap();
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-xdp-test-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("xdp-test.xdp");
    std::fs::write(&path, [buf.clone(), buf].concat()).unwrap();
    let got = load_dir(&dir).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(got.len(), 2);
    let c = &got[1];
    assert!(c.tb && c.check_fence);
    assert_eq!((c.kind, c.n, c.len1, c.len2), (KIND_FLAT, 16, 4, 4));
    assert_eq!((c.gap_open, c.gap_extend, c.x_drop), (11, 1, 38));
    // rows read backwards from base 4: q[3], q[2], q[1], q[0]
    assert_eq!(c.rows, vec![1, 7, 5, 3]);
    // columns read backwards from base 4, zero from k = 4
    assert_eq!(c.cols, vec![10, 8, 6, 4, 0]);
    assert_eq!(c.matrix, m);
    assert_eq!(c.result, Some(res));
    assert_eq!(c.ops, ops.to_vec());
    for a in 1..=4 {
        assert_eq!(c.row_seq().get(a), rows.get(a));
    }
    for k in 0..=6 {
        assert_eq!(c.col_seq().get(k), cols.get(k));
    }
}
