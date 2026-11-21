use super::super::{DiscTemplateType, MaskedInterval};
use super::*;

/// A small deterministic generator (xorshift64*).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// Random BLASTNA queries over few distinct letters (so cells get long chains), with
/// ambiguities and masks.
fn queries(
    rng: &mut Rng,
    n: usize,
    len: usize,
) -> (Vec<Vec<u8>>, Vec<i32>, Vec<Vec<MaskedInterval>>) {
    let mut seqs = Vec::new();
    let mut offsets = Vec::new();
    let mut masks = Vec::new();
    let mut at = 0i32;
    for q in 0..n {
        let l = len / 2 + rng.below(len as u64) as usize;
        let letters = if q % 2 == 0 { 2 } else { 4 };
        let seq: Vec<u8> = (0..l)
            .map(|_| {
                if rng.below(200) == 0 {
                    14 // N, an ambiguity
                } else {
                    rng.below(letters) as u8
                }
            })
            .collect();
        let mut m = Vec::new();
        if q % 3 == 1 && l > 40 {
            let s = rng.below((l - 20) as u64) as usize;
            m.push(MaskedInterval {
                start: s,
                end: s + 15,
            });
        }
        offsets.push(at);
        at += l as i32 + 1;
        seqs.push(seq);
        masks.push(m);
    }
    (seqs, offsets, masks)
}

fn words<'a>(
    seqs: &'a [&'a [u8]],
    offsets: &'a [i32],
    masks: &'a [Vec<MaskedInterval>],
    template_length: usize,
    template_type: DiscTemplateType,
) -> DiscWords<'a> {
    DiscWords {
        queries_blastna: seqs,
        query_offsets: offsets,
        query_masks: masks,
        template_length,
        template_type,
        second_template_type: DiscTemplateType::Contiguous,
        two_templates: false,
    }
}

#[test]
fn the_runs_are_the_chains_with_the_fills_presence_vector_and_counters() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let cases = [
        (DiscTemplateType::T11_16Coding, 16, 22usize, 5usize),
        (DiscTemplateType::T11_18Optimal, 18, 22, 7),
        (DiscTemplateType::T12_21Coding, 21, 24, 8),
    ];
    for (template_type, template_length, log_cells, pv_array_bts) in cases {
        let cells = 1usize << log_cells;
        let pv_len = cells >> pv_array_bts;
        for (n, len) in [(1, 30), (7, 400), (40, 3000)] {
            let (seqs, offsets, masks) = queries(&mut rng, n, len);
            let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
            let w = words(&refs, &offsets, &masks, template_length, template_type);
            let (hashtable, next_pos, fill_pv, fill_helper) =
                chains(&w, cells, pv_len, pv_array_bts);
            let mut pv = vec![0 as PvArrayType; pv_len];
            let mut helper = vec![0u32; cells >> 11];
            let csr = build(&w, cells, &mut pv, pv_array_bts, &mut helper, Vec::new());
            let longest = disc_longest_chain(&fill_helper, &[], false);
            x_assert_runs_are_the_chains(
                &w,
                &csr,
                false,
                &pv,
                &helper,
                pv_array_bts,
                &hashtable,
                &next_pos,
                &fill_pv,
                &fill_helper,
                longest,
            );
            let chained: usize = (0..cells).map(|c| csr.run(c).len()).sum();
            assert_eq!(chained, csr.words());
            assert!(csr.words() > 0 || n == 1);
        }
    }
}

#[test]
fn contains_finds_exactly_the_runs_offsets() {
    let csr = Csr {
        offsets: vec![0, 0, 3, 4, 4],
        x_bits: Vec::new(),
        qs: vec![90, 40, 7, 12],
        x_keep: None,
    };
    for (cell, q, want) in [
        (0, 7, false),
        (1, 90, true),
        (1, 40, true),
        (1, 7, true),
        (1, 41, false),
        (1, 6, false),
        (2, 12, true),
        (2, 7, false),
        (3, 12, false),
    ] {
        assert_eq!(csr.contains(cell, q), want, "{cell} {q}");
    }
}

#[test]
fn the_ring_emits_every_word_once_in_order() {
    let csr = Csr {
        offsets: vec![0, 2, 2, 5],
        x_bits: Vec::new(),
        qs: vec![9, 3, 30, 20, 10],
        x_keep: None,
    };
    for n in [0usize, 1, 15, 16, 17, 31, 32, 33, 100] {
        let mut ring: Ring<'_, (u32, usize)> = Ring::new(&csr);
        let mut out = Vec::new();
        for i in 0..n {
            ring.push((i as u32, i % 3), |w| w.1, &mut |w| out.push(w.0));
        }
        ring.flush(&mut |w| out.push(w.0));
        assert_eq!(out, (0..n as u32).collect::<Vec<_>>());
    }
}

// LOSAT_X_BNCSRPACK: the byte counts of a packed build give each cell its number of words, past
// 255 too, and the presence bits of the non-zero counts and the ranks in cell order, built by
// two owners of 32 blocks each; a single non-zero byte of any value sets its bit only.
#[test]
fn byte_counts_give_the_presence_ranks_and_exact_counts() {
    for k in 0..64 {
        for v in [1u8, 0x7f, 0x80, 0xff] {
            let mut c = [0u8; 64];
            c[k] = v;
            assert_eq!(x_nonzero_bits(&c), 1u64 << k, "byte {k} = {v}");
        }
    }
    assert_eq!(x_nonzero_bits(&[0; 64]), 0);
    assert_eq!(x_nonzero_bits(&[0xff; 64]), u64::MAX);
    let cells = 4096usize;
    let mut cnt = vec![0u8; cells];
    let mut more = HashMap::new();
    let mut want = vec![0u32; cells];
    let mut x = 0x9E37_79B9u32;
    for (cell, w) in want.iter_mut().enumerate() {
        x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        *w = match cell % 9 {
            0 | 4 | 6 => 0,
            1 => 254 + (x >> 16) % 4,
            2 => 1000 + (x >> 16) % 300,
            _ => (x >> 16) % 3,
        };
        for _ in 0..*w {
            // SAFETY: `cell < cells`; one thread.
            unsafe { x_count8(cnt.as_mut_ptr(), &mut more, cell) };
        }
    }
    // LOSAT_X_BNPOPCNT: in each form this CPU runs.
    for pop in XPop::x_forms() {
        let mut bits = vec![0u64; cells >> 5];
        let mut starts = vec![0u32; cells];
        // SAFETY: 64 blocks of 64 cells: `cells` counts, `2 x 64` words, `cells` ranks at most.
        let (present, r) = unsafe {
            let present = x_presence8(bits.as_mut_ptr(), cnt.as_ptr(), 0, 32, pop)
                + x_presence8(bits.as_mut_ptr(), cnt.as_ptr(), 32, 64, pop);
            let r = x_rank_counts8(
                bits.as_mut_ptr(),
                cnt.as_ptr(),
                &more,
                starts.as_mut_ptr(),
                0,
                32,
                0,
            );
            let r = x_rank_counts8(
                bits.as_mut_ptr(),
                cnt.as_ptr(),
                &more,
                starts.as_mut_ptr(),
                32,
                64,
                r,
            );
            (present, r)
        };
        assert_eq!(present, want.iter().filter(|&&w| w > 0).count());
        assert_eq!(r, present);
        let mut rank = 0usize;
        for (cell, &w) in want.iter().enumerate() {
            if w == 0 {
                assert_eq!(x_rank(&bits, cell, pop), None, "cell {cell} {pop:?}");
            } else {
                assert_eq!(x_rank(&bits, cell, pop), Some(rank), "cell {cell} {pop:?}");
                assert_eq!(starts[rank], w, "cell {cell} {pop:?}");
                rank += 1;
            }
        }
    }
}

// LOSAT_X_BNPOPCNT: both forms count the same bits of random and edge words, and rank every
// cell of random presence bits alike (the `popcnt` form where this CPU has it).
#[test]
fn both_counts_give_the_same_ranks() {
    let forms = XPop::x_forms();
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let mut words = vec![
        0u64,
        1,
        u64::MAX,
        1 << 63,
        u64::MAX >> 1,
        0x5555_5555_5555_5555,
    ];
    for _ in 0..10_000 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        words.push(x);
        words.push(x & (x >> 3) & (x >> 5));
    }
    for &w in &words {
        for &pop in &forms {
            assert_eq!(pop.ones(w), w.count_ones() as usize, "{w:#x} {pop:?}");
        }
    }
    let mut bits = Vec::new();
    let mut before = 0u64;
    for &w in &words[..2048] {
        bits.push(w);
        bits.push(before);
        before += u64::from(w.count_ones());
    }
    for cell in 0..(bits.len() / 2) * 64 {
        let want = x_rank(&bits, cell, XPop(false));
        for &pop in &forms {
            assert_eq!(x_rank(&bits, cell, pop), want, "cell {cell} {pop:?}");
        }
    }
}

// LOSAT_X_BNCSRPACK: the packed head of the table above (cells 0 and 2 non-empty, then 61 empty
// cells, then cell 64 with one word) gives the same runs, and its ring emits the words of the
// non-empty cells, in order, and no other; a head of one element per cell has no packed reader.
#[test]
fn a_packed_head_gives_the_same_runs_and_the_ring_skips_empty_cells() {
    let bits = vec![0b101u64, 0, 1, 2];
    let csr = Csr {
        offsets: vec![0, 2, 5, 6],
        x_bits: bits,
        qs: vec![9, 3, 30, 20, 10, 77],
        x_keep: None,
    };
    let packed = csr.x_packed().expect("a packed head");
    assert_eq!(packed.run(0), &[9, 3]);
    assert!(packed.run(1).is_empty() && packed.run(63).is_empty());
    assert_eq!(packed.run(2), &[30, 20, 10]);
    assert_eq!(packed.run(64), &[77]);
    assert_eq!(csr.x_run_any(2), &[30, 20, 10]);
    assert!(csr.x_run_any(1).is_empty());
    let cells = [0usize, 1, 2, 64, 5];
    for n in [0usize, 1, 15, 16, 17, 31, 32, 33, 100] {
        let mut ring: XRingPacked<'_, (u32, usize)> = packed.ring();
        let mut out = Vec::new();
        // LOSAT_X_BNRANKONCE: the element handed over with a word is its run's.
        let mut emit = |w: (u32, usize), e: usize| {
            assert_eq!(packed.run_at(e), packed.run(w.1), "word {}", w.0);
            out.push(w.0);
        };
        for i in 0..n {
            XRing::push(&mut ring, (i as u32, cells[i % 5]), |w| w.1, &mut emit);
        }
        XRing::flush(&mut ring, |w| w.1, &mut emit);
        let want: Vec<u32> = (0..n as u32)
            .filter(|&i| !packed.run(cells[i as usize % 5]).is_empty())
            .collect();
        assert_eq!(out, want);
    }
    let dense = Csr {
        offsets: vec![0, 2, 2, 5],
        x_bits: Vec::new(),
        qs: vec![9, 3, 30, 20, 10],
        x_keep: None,
    };
    assert!(dense.x_packed().is_none());
    assert_eq!(dense.x_run_any(2), dense.run(2));
    // LOSAT_X_BNRANKONCE: the ring of a head of one element per cell emits every word, with
    // its cell as the element.
    for n in [0usize, 1, 31, 32, 33, 100] {
        let mut ring = XHead::ring::<(u32, usize)>(&dense);
        let mut out = Vec::new();
        let mut emit = |w: (u32, usize), e: usize| {
            assert_eq!(e, w.1, "word {}", w.0);
            assert_eq!(XHead::run_at(&dense, e), dense.run(w.1), "word {}", w.0);
            out.push(w.0);
        };
        for i in 0..n {
            XRing::push(&mut ring, (i as u32, i % 3), |w| w.1, &mut emit);
        }
        XRing::flush(&mut ring, |w| w.1, &mut emit);
        assert_eq!(out, (0..n as u32).collect::<Vec<_>>());
    }
}

#[test]
fn the_size_rule_takes_a_quarter_of_the_cells() {
    assert!(!wanted(0, 1 << 22));
    assert!(!wanted((1 << 20) - 1, 1 << 22));
    assert!(wanted(1 << 20, 1 << 22));
    assert!(wanted(usize::MAX, 1 << 22));
}
