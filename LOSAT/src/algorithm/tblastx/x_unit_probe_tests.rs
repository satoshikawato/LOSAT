use super::*;

fn run(p: &mut XStatsProbe, width: usize, batches: &[&[(usize, u32, XHit)]]) -> XStatTotals {
    p.unit_begin(width);
    for batch in batches {
        p.batch(batch.len());
        for &(cell, s, outcome) in batch.iter() {
            p.at(cell, s);
            p.hit(outcome, 40);
        }
        p.batch_end();
    }
    // `unit_end` adds the counts of the cells' last hits and open runs to `p.t`.
    p.unit_end();
    p.t
}

// No NCBI counterpart (statistics). Hits per cell over the unit and per batch, isolated
// hits, and runs, on a hand-made sequence of hits.
#[test]
fn the_statistics_count_cells_isolation_and_runs() {
    use XHit::*;
    let mut p = XStatsProbe::default();
    // Cell 0: hits at 0, 100, 110, 120. Cell 1: one hit. Cell 2: masked, masked, too far,
    // too far, too far, extended.
    let t = run(
        &mut p,
        4,
        &[
            &[
                (0, 0, TooFar),
                (2, 5, Masked),
                (2, 6, Masked),
                (1, 7, TooFar),
            ],
            &[
                (0, 100, TooFar),
                (2, 200, TooFar),
                (2, 300, TooFar),
                (2, 400, TooFar),
            ],
            &[(0, 110, Overlap), (0, 120, Extended), (2, 410, Extended)],
        ],
    );
    assert_eq!(t.hits, 11);
    assert_eq!(t.unit_cells, 3);
    // Unit: cell 0 has 4 hits (4-7), cell 1 one, cell 2 six (4-7).
    assert_eq!(t.unit_cell_hist, [1, 0, 0, 2, 0]);
    assert_eq!(t.unit_hit_hist, [1, 0, 0, 10, 0]);
    // Batches: {0:1, 2:2, 1:1}, {0:1, 2:3}, {0:2, 2:1}.
    assert_eq!(t.batch_cells, 7);
    assert_eq!(t.batch_cell_hist, [4, 2, 1, 0, 0]);
    assert_eq!(t.batch_hit_hist, [4, 4, 3, 0, 0]);
    // Isolated: cell 0 at 0 (100 has 110 near it); cell 1 at 7; cell 2 at 200 and 300 (5
    // and 6 are near each other, and 400 and 410).
    assert_eq!(t.isolated, 4);
    // Runs: cell 0 TooFar x2 then Overlap (one too-far run of 2); cell 1 TooFar x1; cell 2
    // Masked x2, TooFar x3.
    assert_eq!(t.runs[0], [0, 1, 0, 0, 0, 0, 0]);
    assert_eq!(t.runs[1], [1, 1, 1, 0, 0, 0, 0]);
    assert_eq!(t.run_hits, [2, 6]);
}

// No NCBI counterpart (statistics).
#[test]
fn buckets_cover_their_ranges() {
    assert_eq!(
        [1, 2, 3, 4, 7, 8, 1000].map(h_bucket),
        [0, 1, 2, 3, 3, 4, 4]
    );
    assert_eq!(
        [1, 2, 3, 4, 7, 8, 15, 16, 63, 64].map(r_bucket),
        [0, 1, 2, 3, 3, 4, 4, 5, 5, 6]
    );
}
