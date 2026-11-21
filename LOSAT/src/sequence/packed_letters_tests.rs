use super::*;

type Blocks = fn(&[u8], &mut [u8], &mut Rest) -> usize;

fn pack_with(letters: &[u8], blocks: Blocks) -> (Vec<u8>, Vec<OtherRun>, Vec<(u32, u32)>) {
    let mut packed = vec![0u8; letters.len().div_ceil(4)];
    let mut rest = Rest::default();
    let done = blocks(letters, &mut packed, &mut rest);
    for (index, &byte) in letters.iter().enumerate().skip(done) {
        rest.byte(index, byte, &mut packed, index);
    }
    (packed, rest.others, rest.lowercase)
}

fn bytewise(_: &[u8], _: &mut [u8], _: &mut Rest) -> usize {
    0
}

/// A small xorshift, so the cases are the same on every run.
struct Draws(u64);

impl Draws {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

fn cases() -> Vec<Vec<u8>> {
    let mut draws = Draws(0x9E37_79B9_7F4A_7C15);
    let mut cases: Vec<Vec<u8>> = vec![Vec::new(), b"A".to_vec(), b"acgtACGTn".to_vec()];
    let others = b"NnRYKMSWBDHVrykmswbdhvUu-*Xx\x00\xff \x7f";
    for case in 0..3000 {
        let length = (draws.next() % 400) as usize;
        let letters = (0..length)
            .map(|_| match draws.next() % 8 {
                0 if case % 4 != 3 => others[(draws.next() % others.len() as u64) as usize],
                _ if case % 3 == 0 => b"acgt"[(draws.next() & 3) as usize],
                _ if case % 3 == 1 && draws.next().is_multiple_of(64) => {
                    b"acgt"[(draws.next() & 3) as usize]
                }
                _ => b"ACGT"[(draws.next() & 3) as usize],
            })
            .collect();
        cases.push(letters);
    }
    // long records as a reader gives them: runs of N, a lower-case stretch
    let mut long: Vec<u8> = (0..5000)
        .map(|_| b"ACGT"[(draws.next() & 3) as usize])
        .collect();
    long[1000..1037].fill(b'N');
    long[2001..2300].make_ascii_lowercase();
    cases.push(long);
    cases
}

#[test]
fn the_letters_are_made_again_exactly() {
    let mut packed_records = 0;
    for letters in cases() {
        let mut out = b"left over".to_vec();
        if let Some(packed) = PackedLetters::pack(&letters) {
            assert_eq!(packed.len(), letters.len());
            packed.unpack_into(&mut out);
            assert_eq!(out, letters);
            packed_records += 1;
        }
    }
    assert!(packed_records > 400, "{packed_records} records packed");
}

#[test]
fn a_store_in_a_buffer_gives_the_letters_back() {
    let mut buffer = b"head".to_vec();
    let mut stored = Vec::new();
    let mut cases = cases();
    // Long runs of other bytes and long lower-case stretches.
    cases.push([b"ACGT".repeat(40), vec![b'N'; 300], b"acgt".repeat(50)].concat());
    cases.push(b"acgtnnnnACGTNNNN".repeat(30));
    for (index, letters) in cases.iter().enumerate() {
        let masks: &[(u32, u32)] = if index % 5 == 1 && letters.len() > 12 {
            &[(2, 5), (9, 11)]
        } else {
            &[]
        };
        let upper = letters.to_ascii_uppercase();
        let source = if masks.is_empty() { letters } else { &upper };
        if let Some(store) = PackedLetters::pack_masked(source, masks) {
            let start = buffer.len();
            let at = store.append_to(&mut buffer);
            assert_eq!(buffer.len() - start, at.bytes(), "case {index}");
            assert_eq!(store.stored_bytes(), at.bytes(), "case {index}");
            assert_eq!(at.len(), store.len());
            let mut expected = Vec::new();
            store.unpack_into(&mut expected);
            stored.push((start, at, expected));
        }
    }
    assert!(stored.len() > 400, "{} stores", stored.len());
    let mut out = b"old".to_vec();
    for (start, at, expected) in &stored {
        unpack_from(&buffer[*start..start + at.bytes()], *at, &mut out);
        assert_eq!(&out, expected);
    }
    assert_eq!(&buffer[..4], b"head");
}

#[test]
fn every_form_packs_as_byte_by_byte() {
    for letters in cases() {
        let expected = pack_with(&letters, bytewise);
        assert_eq!(
            pack_with(&letters, pack_blocks_swar),
            expected,
            "{letters:?}"
        );
        #[cfg(target_arch = "x86_64")]
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: `pack_blocks_avx2` needs only AVX2, detected by this `if`; `avx2` is
            // called only by `pack_with` in the next statement, inside the same `if`.
            let avx2: Blocks =
                |letters, packed, rest| unsafe { pack_blocks_avx2(letters, packed, rest) };
            assert_eq!(pack_with(&letters, avx2), expected, "{letters:?}");
        }
    }
}

#[test]
fn every_form_unpacks_as_the_table() {
    for letters in cases() {
        let Some(store) = PackedLetters::pack(&letters) else {
            continue;
        };
        let mut table = vec![std::mem::MaybeUninit::new(0u8); letters.len()];
        let mut done = 0;
        for (dst, &byte) in table.chunks_exact_mut(4).zip(&store.packed[..]) {
            for (slot, &letter) in dst.iter_mut().zip(&LETTERS[usize::from(byte)]) {
                slot.write(letter);
            }
            done += 4;
        }
        // SAFETY: every slot of `table` was made by `MaybeUninit::new(0u8)` and then only
        // overwritten by `slot.write`, so each slot read here is initialized.
        let table: Vec<u8> = table[..done]
            .iter()
            .map(|slot| unsafe { slot.assume_init() })
            .collect();
        let mut out = vec![std::mem::MaybeUninit::new(0u8); letters.len()];
        let written = unpack_blocks(&store.packed, &mut out);
        assert!(written <= done && done - written < 32);
        // SAFETY: every slot of `out` was made by `MaybeUninit::new(0u8)`, and `unpack_blocks`
        // only writes letters into it, so each slot read here is initialized.
        let out: Vec<u8> = out[..written]
            .iter()
            .map(|slot| unsafe { slot.assume_init() })
            .collect();
        assert_eq!(out, table[..written]);
    }
}

#[test]
fn a_mask_lower_cases_the_reader_letters() {
    let mut letters = b"ACGTNNACGTACACGTACGTACGTACGTACGTACGTACGTACGTACGT".repeat(8);
    letters.extend_from_slice(b"RYACG");
    let masked = [(2u32, 9u32), (40, 47), (300, 386)];
    let mut expected = letters.clone();
    for &(from, to) in &masked {
        expected[from as usize..=to as usize].make_ascii_lowercase();
    }
    let store = PackedLetters::pack_masked(&letters, &masked).expect("packs");
    let mut out = Vec::new();
    store.unpack_into(&mut out);
    assert_eq!(out, expected);
}

#[test]
fn a_record_that_packs_to_no_fewer_bytes_keeps_its_letters() {
    let letters: Vec<u8> = (0..400).map(|index| b"NR"[index % 2]).collect();
    assert_eq!(PackedLetters::pack(&letters), None);
    let store = PackedLetters::pack(&[b'N'; 400]).expect("one run");
    assert_eq!(store.heap_bytes(), 100 + 12);
    #[cfg(target_pointer_width = "64")]
    assert_eq!(std::mem::size_of::<PackedLetters>(), 32);
    // The sizes below are the 64-bit ones (`HELD_BESIDE` and the chunks follow the pointer
    // width, so a 32-bit target packs records of fewer letters).
    #[cfg(target_pointer_width = "64")]
    {
        // a 32-byte chunk of letters against 40 bytes beside a store and its 32-byte chunk
        assert_eq!(PackedLetters::pack(&[b'A'; 22]), None);
        assert_eq!(PackedLetters::pack(&[b'A'; 50]), None);
        // 80 bytes against 40 + 32: from about 60 letters a record packs
        assert!(PackedLetters::pack(&[b'A'; 64]).is_some());
        assert!(PackedLetters::pack(&[b'A'; 150]).is_some());
    }
}
