use super::*;

#[test]
fn test_simple_repeat() {
    let masker = DustMasker::with_defaults();

    // Simple repeat sequence should be masked
    let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let intervals = masker.mask_sequence(seq);

    // Should have at least one masked interval
    assert!(!intervals.is_empty(), "Poly-A sequence should be masked");
}

#[test]
fn test_complex_sequence() {
    let masker = DustMasker::with_defaults();

    // Use a de Bruijn sequence B(4,3) over A/C/G/T: each 3-mer appears exactly once.
    // This is a more appropriate \"high complexity\" control than a short periodic repeat.
    fn debruijn_acgt_k3() -> Vec<u8> {
        fn db(t: usize, p: usize, k: usize, n: usize, a: &mut [usize], seq: &mut Vec<usize>) {
            if t > n {
                if n.is_multiple_of(p) {
                    for i in 1..=p {
                        seq.push(a[i]);
                    }
                }
            } else {
                a[t] = a[t - p];
                db(t + 1, p, k, n, a, seq);
                for j in (a[t - p] + 1)..k {
                    a[t] = j;
                    db(t + 1, t, k, n, a, seq);
                }
            }
        }

        let alphabet: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let k = alphabet.len();
        let n = 3usize;
        let mut a = vec![0usize; k * n + 1];
        let mut idx_seq: Vec<usize> = Vec::new();
        db(1, 1, k, n, &mut a, &mut idx_seq);

        let mut out: Vec<u8> = idx_seq.into_iter().map(|i| alphabet[i]).collect();
        // Linearize the cyclic de Bruijn sequence by appending the first n-1 symbols.
        let prefix: Vec<u8> = out[..(n - 1)].to_vec();
        out.extend_from_slice(&prefix);
        out
    }

    let seq = debruijn_acgt_k3();
    let intervals = masker.mask_sequence(&seq);

    // Complex sequence should have few or no masked regions
    let total_masked: usize = intervals.iter().map(|i| i.end - i.start).sum();
    assert!(
        total_masked < seq.len() / 2,
        "Complex sequence should not be heavily masked"
    );
}

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:214-216,229-233
/// The test compares the ring-buffer masker with the port of this function on random
/// sequences, levels, windows and ranges.
/// The ring-buffer window must give the reference intervals on random
/// sequences of every kind the masker distinguishes: plain, biased,
/// tandem repeats of every short period, homopolymer runs and Ns.
#[test]
fn x_ring_buffer_dust_matches_reference_on_random_sequences() {
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3000);
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let bases = b"ACGT";
    let mut masked_cases = 0usize;
    for case in 0..cases {
        let len = 1 + rng.below(1500) as usize;
        let mut seq = Vec::with_capacity(len);
        while seq.len() < len {
            let piece = 1 + rng.below(200) as usize;
            match rng.below(7) {
                0 | 1 => {
                    for _ in 0..piece {
                        seq.push(bases[rng.below(4) as usize]);
                    }
                }
                2 => {
                    // biased composition
                    let major = bases[rng.below(4) as usize];
                    for _ in 0..piece {
                        seq.push(if rng.below(5) == 0 {
                            bases[rng.below(4) as usize]
                        } else {
                            major
                        });
                    }
                }
                3 => {
                    // tandem repeat with a few substitutions
                    let period = 1 + rng.below(7) as usize;
                    let unit: Vec<u8> = (0..period).map(|_| bases[rng.below(4) as usize]).collect();
                    for i in 0..piece {
                        seq.push(if rng.below(25) == 0 {
                            bases[rng.below(4) as usize]
                        } else {
                            unit[i % period]
                        });
                    }
                }
                4 => {
                    let base = bases[rng.below(4) as usize];
                    seq.extend(std::iter::repeat_n(base, piece));
                }
                5 => {
                    let n = 1 + rng.below(12) as usize;
                    seq.extend(std::iter::repeat_n(b'N', n));
                }
                _ => {
                    // other IUPAC letters and lower case
                    for _ in 0..piece.min(20) {
                        seq.push(b"acgtRYKMSWnBDHV"[rng.below(15) as usize]);
                    }
                }
            }
        }
        seq.truncate(len);
        let (level, window, linker) = match case % 4 {
            0 => (20, 64, 1),
            1 => (
                10 + rng.below(40) as u32,
                8 + rng.below(57) as usize,
                1 + rng.below(32) as usize,
            ),
            2 => (20, 64, 1),
            _ => (2 + rng.below(63) as u32, 64, 1),
        };
        let masker = DustMasker::new(level, window, linker);
        let start = if case % 5 == 0 {
            rng.below(len as u64) as usize
        } else {
            0
        };
        // `stop` is the last position, `start` at most `stop` (as `mask_range` makes them).
        let stop = if case % 7 == 0 {
            rng.below(len as u64) as usize
        } else {
            len - 1
        };
        let start = start.min(stop);
        let reference = masker.mask_range_reference(&seq, start, stop, &mut NcbiRandom::default());
        let fast = masker.x_mask_range(&seq, start, stop, &mut NcbiRandom::default());
        assert_eq!(
            fast, reference,
            "case {case}: level={level} window={window} linker={linker} start={start} stop={stop}"
        );
        // `LOSAT_X_DUSTREG` (`x_dust_regs.rs`): the same intervals as the reference.
        let regs = masker.x_mask_range_regs(&seq, start, stop, &mut NcbiRandom::default());
        assert_eq!(
            regs, reference,
            "regs case {case}: level={level} window={window} linker={linker} start={start} stop={stop}"
        );
        if !reference.is_empty() {
            masked_cases += 1;
        }
    }
    // the generator must actually exercise the masking paths
    assert!(
        masked_cases * 2 > cases,
        "only {masked_cases} of {cases} cases masked anything"
    );
}

#[test]
fn test_short_sequence() {
    let masker = DustMasker::with_defaults();

    // Very short sequences should return empty
    let seq = b"AC";
    let intervals = masker.mask_sequence(seq);
    assert!(intervals.is_empty());
}

#[test]
fn test_params_validation() {
    // Test parameter validation
    let params = DustParams::new(1, 5, 0);
    assert_eq!(params.level, 20); // Should default to 20 (out of range)
    assert_eq!(params.window, 64); // Should default to 64 (out of range)
    assert_eq!(params.linker, 1); // Should default to 1 (out of range)

    let params = DustParams::new(30, 32, 16);
    assert_eq!(params.level, 30);
    assert_eq!(params.window, 32);
    assert_eq!(params.linker, 16);
}
