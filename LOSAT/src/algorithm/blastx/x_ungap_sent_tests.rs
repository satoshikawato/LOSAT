use super::*;

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

/// A residue code: mostly 1..=24, sometimes the stop codon 25, a code 26..=31, or (with
/// `zero`) an inner 0.
fn code(rng: &mut Rng, zero: bool) -> u8 {
    match rng.below(40) {
        0 => 25,
        1 => 26 + rng.below(6) as u8,
        2 if zero => 0,
        _ => 1 + rng.below(24) as u8,
    }
}

/// A six-frame query buffer as `query_setup.rs` lays it out: `[0, frame, 0, frame, .., 0]`,
/// frames of 1 letter and more, some copied from `subject` with substitutions so that walks
/// run long.
fn query_buffer(rng: &mut Rng, subject: &[u8]) -> Vec<u8> {
    let mut buffer = vec![0u8];
    for _ in 0..6 {
        let most = if rng.below(4) == 0 { 3 } else { 90 };
        let len = 1 + rng.below(most) as usize;
        if rng.below(2) == 0 && subject.len() > 2 {
            let from = 1 + rng.below((subject.len() - 2) as u64) as usize;
            for i in 0..len {
                let c = subject[1 + (from + i) % (subject.len() - 2)];
                let c = if rng.below(6) == 0 {
                    code(rng, false)
                } else {
                    c
                };
                buffer.push(c);
            }
        } else {
            for _ in 0..len {
                buffer.push(code(rng, false));
            }
        }
        buffer.push(0);
    }
    buffer
}

/// The path gives NCBI's extension (`extend_two_hit`, the reference)
/// on the reference's sequences (the buffers without their first byte) for every offset of a
/// word, the ends of both sequences and the frame borders included, subjects of 0 to 3
/// letters, inner sentinels, codes up to 31, and X-drops on and off the path (above 127: the
/// reference).
#[test]
fn two_hit_is_the_reference_on_blastx_buffers() {
    two_hit_is_the_reference_with::<3>(ScoringMatrix::Blosum62);
}

/// [`two_hit_is_the_reference_on_blastx_buffers`] for the word sizes 2 and 4 (item SC9b-1):
/// the words of 2 and 4 letters at the ends of both sequences, over inner and last sentinels,
/// BLOSUM62 and PAM30 (the largest scores of the matrices).
#[test]
fn two_hit_of_word_sizes_2_and_4_is_the_reference() {
    for matrix in [ScoringMatrix::Blosum62, ScoringMatrix::Pam30] {
        two_hit_is_the_reference_with::<2>(matrix);
        two_hit_is_the_reference_with::<4>(matrix);
    }
}

/// [`two_hit_is_the_reference_on_blastx_buffers`] with every other matrix (item SC6b): all of
/// them have keys (IDENTITY too, which ungapped searches take), and a query batch of theirs is
/// on the path where the switch is on.
#[test]
fn two_hit_of_other_matrices_is_the_reference() {
    for matrix in ScoringMatrix::ALL {
        if matrix == ScoringMatrix::Blosum62 {
            continue;
        }
        assert!(keys_of(matrix).is_some(), "{matrix:?} has keys");
        for word in [2, 3, 4] {
            assert_eq!(
                query_fits_of(matrix, &[0, 3, 0], word, 40),
                on(),
                "{matrix:?} word {word}"
            );
        }
        for word in [1, 5, 6, 7] {
            assert!(
                !query_fits_of(matrix, &[0, 3, 0], word, 40),
                "{matrix:?} word {word}"
            );
        }
        assert!(
            !query_fits_of(matrix, &[0, 3, 0], 3, 0),
            "{matrix:?} one-hit"
        );
        two_hit_is_the_reference_with::<3>(matrix);
    }
}

/// The path with the word size `W` against `extend_two_hit` with the same word size.
fn two_hit_is_the_reference_with<const W: usize>(matrix: ScoringMatrix) {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15 ^ matrix as u64 ^ ((W as u64) << 40));
    let mut compared = 0u64;
    for round in 0..400 {
        let slen = rng.below(if round % 8 == 0 { 4 } else { 120 }) as usize;
        let mut subject = vec![0u8];
        subject.extend((0..slen).map(|_| code(&mut rng, true)));
        subject.push(0);
        let query = query_buffer(&mut rng, &subject);
        // The walks of BLOSUM62's type, or of another matrix's.
        let blosum62 = (matrix == ScoringMatrix::Blosum62)
            .then(|| XSent::for_test(super::super::Blosum62Bx, &query, &subject, false));
        let other = (matrix != ScoringMatrix::Blosum62)
            .then(|| XSent::for_test(super::super::MatrixBx::of(matrix), &query, &subject, false));
        let (q_ref, s_ref) = (&query[1..], &subject[1..]);
        for _ in 0..200 {
            let q = rng.below(q_ref.len() as u64 + 2) as usize;
            let s = rng.below(s_ref.len() as u64 + 2) as usize;
            let s_left = rng.below(s as u64 + 4) as i32;
            let dropoff = [1, 7, 16, 40, 127, 128, 300][rng.below(7) as usize];
            let want = extend_two_hit(matrix, q_ref, s_ref, s_left, s, q, dropoff, W);
            let got = match (&blosum62, &other) {
                (Some(sent), _) => sent.two_hit::<false, W>(s_left, s, q, dropoff),
                (_, Some(sent)) => sent.two_hit::<false, W>(s_left, s, q, dropoff),
                _ => unreachable!(),
            };
            assert_eq!(
                got, want,
                "{matrix:?} W={W} round {round}: s_left={s_left} s={s} q={q} dropoff={dropoff}"
            );
            compared += u64::from(want.is_some());
        }
    }
    assert!(compared > 10_000, "only {compared} extensions");
}

/// The batch check: sentinels at both ends and codes below 32 (`query_fits` adds the switch,
/// the word size and the window).
#[test]
fn query_on_path_checks_the_buffer() {
    assert!(word_size_on_path(2) && word_size_on_path(3) && word_size_on_path(4));
    assert!(!word_size_on_path(1) && !word_size_on_path(5) && !word_size_on_path(7));
    assert!(query_on_path(&[0, 3, 4, 0, 5, 0]));
    assert!(query_on_path(&[0, 0]));
    assert!(!query_on_path(&[0]));
    assert!(!query_on_path(&[]));
    assert!(!query_on_path(&[1, 3, 0]));
    assert!(!query_on_path(&[0, 3, 1]));
    assert!(!query_on_path(&[0, 32, 0]));
    assert!(query_on_path(&[0, 31, 0]));
}
