use super::blast_get_start_for_gapped_alignment_nucl;

// `x_same_matrix_fold` (wasm32's `x_same_matrix`) against the slice comparison: equal matrices,
// and a matrix that differs from the cached one in one cell by one bit (every cell, every bit),
// or in its sign.
#[test]
fn x_same_matrix_is_the_slice_comparison() {
    use super::{x_same_matrix, x_same_matrix_fold, BlastnaMatrix};
    let mut cached: BlastnaMatrix = [0; 256];
    for (k, v) in cached.iter_mut().enumerate() {
        *v = (k as i32 % 7) - 3;
    }
    cached[15] = i32::MIN / 2;
    let same = cached;
    assert!(x_same_matrix_fold(&cached, &same));
    assert!(x_same_matrix(&cached, &same));
    assert_eq!(x_same_matrix_fold(&cached, &same), cached[..] == same[..]);
    for cell in 0..256 {
        for bit in 0..32 {
            let mut other = cached;
            other[cell] ^= 1i32 << bit;
            assert!(
                !x_same_matrix_fold(&cached, &other),
                "cell {cell} bit {bit}"
            );
            assert!(!x_same_matrix(&cached, &other), "cell {cell} bit {bit}");
            assert_eq!(x_same_matrix_fold(&cached, &other), cached[..] == other[..]);
        }
        let mut other = cached;
        other[cell] = other[cell].wrapping_neg().wrapping_sub(1);
        assert_eq!(x_same_matrix_fold(&cached, &other), cached[..] == other[..]);
    }
}

// A final X-drop far past the sequences (`-xdrop_gap_final 1e9`, raw 1.1e9 with 2/-3 5/2)
// gives the report of a saturated one (1e4 bits, NCBI 2.17.0 line RG X1b probe `b2.A.*`)
// and holds only the cells that NCBI touches: `N + 1` DP cells (or the first 1000) and rows of
// at most `N + 2` cells, not `x_dropoff / gap_extend` of each.
#[test]
fn a_large_final_x_drop_holds_only_the_sequences_cells() {
    let mut state = 12345u32;
    let mut base = || {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
        ((state >> 16) & 3) as u8
    };
    let q: Vec<u8> = (0..300).map(|_| base()).collect();
    let mut s = q.clone();
    for i in (7..300).step_by(11) {
        s[i] = (s[i] + 1) & 3;
    }
    s.drain(150..153);
    let n = s.len();
    let matrix = super::build_blastna_matrix(2, -3);
    let run = |x_drop: i32| {
        let mut scratch = super::GapAlignScratch::new();
        let out = super::extend_gapped_one_direction_with_traceback_with_scratch(
            &q,
            &s,
            2,
            -3,
            &matrix,
            5,
            2,
            x_drop,
            &mut scratch,
        );
        (out, scratch)
    };
    let (saturated, _) = run(11_088);
    let (large, scratch) = run(1_100_000_000);
    assert_eq!(large, saturated);
    assert_eq!(scratch.dp_mem.len(), 1000);
    // No NCBI counterpart: the vector kernels run this X-drop as the clamp's `R`
    // (`xdrop_simd::x_drop_for_lanes`, `LOSAT_X_DP32`) and leave the scalar port's count as it
    // was; the scalar port itself, which answers without them (other targets,
    // `LOSAT_X_DP32=0`) and for the calls they decline, keeps NCBI's count with only the
    // sequences' cells (below).
    if crate::utils::xdrop_simd::x_kernels_run() && crate::utils::xdrop_simd::x_dp32() {
        assert_eq!(scratch.dp_mem_alloc, 1000);
    }
    let run_scalar = |x_drop: i32| {
        let mut scratch = super::GapAlignScratch::new();
        let out = super::extend_gapped_one_direction_with_traceback_scalar(
            &q,
            &s,
            2,
            -3,
            &matrix,
            5,
            2,
            x_drop,
            &mut scratch,
        );
        (out, scratch)
    };
    let (large, scratch) = run_scalar(1_100_000_000);
    assert_eq!(large, saturated);
    assert_eq!(scratch.dp_mem.len(), 1000);
    assert!(scratch.dp_mem_alloc > 500_000_000);
    assert!(scratch.trace_rows.iter().all(|row| row.len() <= n + 2));
    assert!(scratch
        .trace_rows
        .iter()
        .all(|row| row.capacity() <= 3 * (n + 2)));
}

fn blastna(bases: &str) -> Vec<u8> {
    bases
        .bytes()
        .map(|base| match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
        })
        .collect()
}

// The backward identity run reaches the first subject base before the first
// query base: NCBI stops on the subject's leading sentinel.
#[test]
fn gapped_start_stops_at_the_subject_start() {
    let query = blastna("GGGGGACGT");
    let subject = blastna("ACGTCCCC");
    let query = crate::algorithm::blastn::alignment::ContextQuery {
        buf: &query,
        start: 0,
        len: query.len() as i32,
    };
    let start = blast_get_start_for_gapped_alignment_nucl(query, &subject, 5, 9, 0, 4, 5, 0);
    assert_eq!(start, (7, 2));
}

// The forward identity run reaches the subject end before the query end:
// NCBI stops on the subject's trailing sentinel.
#[test]
fn gapped_start_stops_at_the_subject_end() {
    let query = blastna("ACGTACGTAA");
    let subject = blastna("ACGTACGT");
    let query = crate::algorithm::blastn::alignment::ContextQuery {
        buf: &query,
        start: 0,
        len: query.len() as i32,
    };
    let start = blast_get_start_for_gapped_alignment_nucl(query, &subject, 0, 10, 0, 8, 0, 0);
    assert_eq!(start, (4, 4));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2986-3003
// Both sides of `LOSAT_X_DP32` in one process: with X-drops past the 16-bit lanes (50,000, 1e6
// and 2^31 - 300) and scores large enough for the 32-bit lanes on short sequences (reward
// 1,000, penalty -1,500, gap costs 2,500 / 1,000), the dispatched extensions (the packed
// score-only `s_BlastAlignPackedNucl`, left and right of the seed as
// `extend_gapped_heuristic_with_scratch` calls them, and the traceback `ALIGN_EX` forward:
// the 32-bit lanes) give what the scalar ports give (`blast_align_packed_nucl_scalar`,
// `extend_gapped_one_direction_with_traceback_scalar`: what `LOSAT_X_DP32=0` runs).
#[test]
fn wide_lanes_give_the_scalar_ports_results() {
    let matrix = super::build_blastna_matrix(1000, -1500);
    let (gap_open, gap_extend) = (2500, 1000);
    let mut state = 0x9E37_79B9u32;
    let mut next = move |n: u32| {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        (state >> 8) % n
    };
    let runs = crate::utils::xdrop_simd::X32_RUNS.with(|n| n.get());
    for case in 0..12 {
        let len = 80 + next(400) as usize;
        let q: Vec<u8> = (0..len).map(|_| next(4) as u8).collect();
        let mut s = Vec::with_capacity(len + 8);
        for (i, &c) in q.iter().enumerate() {
            match next(50) {
                0 => s.extend([c, next(4) as u8]),
                1 if i > 0 => {}
                2..=5 => s.push(next(4) as u8),
                _ => s.push(c),
            }
        }
        // NCBI2NA: four bases a byte, the first in the high bits
        let mut packed = vec![0u8; s.len() / 4 + 1];
        for (k, &b) in s.iter().enumerate() {
            packed[k / 4] |= b << (2 * (3 - k % 4));
        }
        let (qs, ss) = (len / 2, (len / 2).min(s.len() - 8));
        let adjustment = super::COMPRESSION_RATIO - (ss % super::COMPRESSION_RATIO);
        let (q_length, s_length) = (qs + adjustment, ss + adjustment);
        for x_drop in [50_000, 1_000_000, i32::MAX - 300] {
            let packed_calls = [
                (q_length, s_length, 0, true, 0),
                (
                    q.len() - q_length,
                    s.len() - s_length,
                    q_length as isize - 1,
                    false,
                    ((s_length + 3) / super::COMPRESSION_RATIO) as isize - 1,
                ),
            ];
            for (query_len, subject_len, query_offset, reverse, byte_offset) in packed_calls {
                let (mut a, mut b) = (super::GapAlignScratch::new(), super::GapAlignScratch::new());
                let got = super::blast_align_packed_nucl_with_scratch(
                    &q,
                    &packed,
                    query_len,
                    subject_len,
                    query_offset,
                    &matrix,
                    gap_open,
                    gap_extend,
                    x_drop,
                    reverse,
                    byte_offset,
                    &mut a,
                );
                let want = super::blast_align_packed_nucl_scalar(
                    &q,
                    &packed,
                    query_len,
                    subject_len,
                    query_offset,
                    &matrix,
                    gap_open,
                    gap_extend,
                    x_drop,
                    reverse,
                    byte_offset,
                    &mut b,
                );
                assert_eq!(
                    got, want,
                    "case {case} x {x_drop} reverse {reverse}: score only"
                );
            }
            let (mut a, mut b) = (super::GapAlignScratch::new(), super::GapAlignScratch::new());
            let got = super::extend_gapped_one_direction_with_traceback_with_scratch(
                &q[qs..],
                &s[ss..],
                1000,
                -1500,
                &matrix,
                gap_open,
                gap_extend,
                x_drop,
                &mut a,
            );
            let want = super::extend_gapped_one_direction_with_traceback_scalar(
                &q[qs..],
                &s[ss..],
                1000,
                -1500,
                &matrix,
                gap_open,
                gap_extend,
                x_drop,
                &mut b,
            );
            assert_eq!(got, want, "case {case} x {x_drop}: traceback");
        }
    }
    if crate::utils::xdrop_simd::x_kernels_run() && crate::utils::xdrop_simd::x_dp32() {
        let runs = crate::utils::xdrop_simd::X32_RUNS.with(|n| n.get()) - runs;
        assert!(runs > 80, "the 32-bit lanes ran {runs} times");
    }
}
