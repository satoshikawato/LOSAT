//! Tests of the BitMagic decoder. The streams named `w..`, `s..`, `m..` are the Q-BM probes of
//! line RG; their expected sizes and bits are what NCBI's BitMagic 8.0.1 decoded (a harness
//! built from NCBI's headers, outside this repository). `differential` is the driver of the
//! differential test against that harness on generated streams: it does nothing unless
//! `BM_DIFF_IN` names a case file.

use super::*;

fn hex(text: &str) -> Vec<u8> {
    let text: String = text.chars().filter(|c| !c.is_whitespace()).collect();
    (0..text.len() / 2)
        .map(|i| u8::from_str_radix(&text[2 * i..2 * i + 2], 16).unwrap())
        .collect()
}

fn ok(size: u32, data: &[u8]) -> BitVector {
    match decode(size, data) {
        Ok(vector) => vector,
        Err(fail) => panic!("{fail:?}"),
    }
}

fn fail(size: u32, data: &[u8]) -> Fail {
    match decode(size, data) {
        Ok(vector) => panic!("decoded {:?}", vector.ranges()),
        Err(fail) => fail,
    }
}

/// A stream with neither byte order nor GAP levels in the header (`BM_HM_NO_BO | BM_HM_NO_GAPL`).
fn stream(tokens: &[u8]) -> Vec<u8> {
    let mut data = vec![0x18];
    data.extend_from_slice(tokens);
    data
}

/// `bm::bit_out` (encoding.h:1301-1360, 1532-1560, 1644-1680): bits lowest first, flushed as
/// little-endian 32-bit words.
#[derive(Default)]
struct BitOut {
    words: Vec<u32>,
    acc: u64,
    used: u32,
}

impl BitOut {
    fn bits(&mut self, value: u32, count: u32) {
        let value = u64::from(value) & ((1u64 << count) - 1);
        self.acc |= value << self.used;
        self.used += count;
        while self.used >= 32 {
            self.words.push(self.acc as u32);
            self.acc >>= 32;
            self.used -= 32;
        }
    }

    fn gamma(&mut self, value: u32) {
        assert!(value > 0);
        let logv = 31 - value.leading_zeros();
        self.bits(0, logv);
        self.bits(1, 1);
        if logv > 0 {
            self.bits(value & ((1 << logv) - 1), logv);
        }
    }

    fn bic(&mut self, arr: &[u32], lo: u32, hi: u32) {
        let mut arr = arr;
        let mut lo = lo;
        while !arr.is_empty() {
            let sz = arr.len() as u32;
            let mid = sz >> 1;
            let val = arr[mid as usize];
            let r = hi.wrapping_sub(lo).wrapping_sub(sz).wrapping_add(1);
            if r != 0 {
                let x = val.wrapping_sub(lo).wrapping_sub(mid);
                let n = r.wrapping_add(1);
                let mut logv = 31 - n.leading_zeros();
                let c = ((1u64 << (logv + 1)) - u64::from(n)) as u32;
                let half_c = i64::from(c >> 1);
                let half_r = i64::from(r >> 1);
                let lo1 = half_r - half_c - i64::from(n & 1);
                let hi1 = half_r + half_c + 1;
                logv += u32::from(i64::from(x) <= lo1 || i64::from(x) >= hi1);
                self.bits(x, logv);
            }
            self.bic(&arr[..mid as usize], lo, val.wrapping_sub(1));
            arr = &arr[mid as usize + 1..];
            lo = val.wrapping_add(1);
        }
    }

    fn bytes(mut self) -> Vec<u8> {
        if self.used > 0 {
            self.words.push(self.acc as u32);
        }
        self.words.iter().flat_map(|w| w.to_le_bytes()).collect()
    }
}

fn reader(data: &[u8]) -> Bytes<'_> {
    Bytes {
        data,
        pos: 0,
        beq: false,
    }
}

#[test]
fn the_header_is_read_in_ncbis_order() {
    // `&GetData()[0]` of an empty OCTET STRING: NCBI reads its first byte from nowhere.
    assert_eq!(fail(0, &[]), Fail::Refuse(Rule::ReadPastEnd));
    // The byte order is read before the 64-bit flag is looked at.
    assert_eq!(fail(0, &[0x20]), Fail::Refuse(Rule::ReadPastEnd));
    assert_eq!(fail(0, &[0x20, 0x01]), Fail::Invalid);
    assert_eq!(fail(0, &[0x38]), Fail::Invalid);
    // GAP levels, then the size, which only grows the vector.
    let vector = ok(100, &[0x08, 1, 0, 2, 0, 3, 0, 4, 0, 0x00]);
    assert_eq!((vector.size(), vector.count()), (100, 0));
    let vector = ok(100, &[0x1A, 50, 0, 0, 0, 0x00]);
    assert_eq!(vector.size(), 100);
    let vector = ok(100, &[0x1A, 0, 1, 0, 0, 0x00]);
    assert_eq!(vector.size(), 256);
    assert_eq!(fail(0, &[0x1A, 0, 1, 0]), Fail::Refuse(Rule::ReadPastEnd));
    // A byte order other than 1 is NCBI's big-endian decoder.
    let vector = ok(0, &[0x12, 0x00, 0, 0, 1, 0, 0x00]);
    assert_eq!(vector.size(), 256);
}

#[test]
fn an_id_list_sets_its_bits_and_grows_the_size() {
    let data = [
        0x0C, 2, 0, 0, 0, 5, 0, 0, 0, 0xFF, 0xFF, 0xFF, 0xFF, 0xAA, 0xBB,
    ];
    let (vector, end) = decode_to_end(3, &data).unwrap();
    assert_eq!(vector.size(), u32::MAX);
    assert_eq!(vector.ranges(), [(5, 5), (u32::MAX, u32::MAX)]);
    assert_eq!(end, 12);
    let vector = ok(3, &[0x0E, 0, 1, 0, 0, 1, 0, 0, 0, 7, 0, 0, 0]);
    assert_eq!((vector.size(), vector.ranges()), (256, vec![(7, 7)]));
    assert_eq!(
        fail(3, &[0x0C, 2, 0, 0, 0, 5, 0, 0, 0]),
        Fail::Refuse(Rule::ReadPastEnd)
    );
}

#[test]
fn block_runs_and_full_blocks() {
    // 1one at block 0, three zero blocks, 8one with 2 blocks (blocks 4 and 5), end
    let vector = ok(10, &stream(&[2, 3, 3, 4, 2, 0]));
    assert_eq!(vector.size(), 6 << 16);
    assert_eq!(vector.ranges(), [(0, 0xFFFF), (4 << 16, (6 << 16) - 1)]);
    assert_eq!(vector.count(), 3 << 16);
    assert_eq!(vector.count_to(5 << 16).unwrap(), (2 << 16) + 1);
    // aone from block 65534 leaves the size
    let vector = ok(10, &stream(&[5, 0xFE, 0xFF, 10]));
    assert_eq!((vector.size(), vector.count()), (10, 2 << 16));
    assert_eq!(vector.count_to(u32::MAX).unwrap(), 2 << 16);
    assert!(vector.get_bit(u32::MAX).unwrap());
    // 8one with 0 blocks: `set_range(from, from - 1)` sets both bits, the size to `from + 1`
    let vector = ok(0, &stream(&[0x81, 4, 0, 0]));
    assert_eq!(vector.ranges(), [(0xFFFF, 0x10000)]);
    assert_eq!(vector.size(), 0x10001);
    // ... and at block 0, `set_range(0, 2^32 - 1)`
    let vector = ok(0, &stream(&[4, 0, 0]));
    assert_eq!((vector.size(), vector.count()), (u32::MAX, 1 << 32));
    // the short zero run passing the last block ends the loop
    let vector = ok(0, &stream(&[7, 0xFF, 0xFF, 0, 0, 0x80 | 2, 2]));
    assert_eq!(vector.count(), 0);
}

#[test]
fn unknown_tokens_throw_and_xor_tokens_are_refused() {
    for token in [12u8, 13, 25, 46, 70, 100, 127] {
        assert_eq!(fail(0, &stream(&[token])), Fail::Invalid, "{token}");
    }
    assert_eq!(
        fail(0, &stream(&[26, 1, 2, 3, 4, 5, 6, 7])),
        Fail::Refuse(Rule::ReadPastEnd)
    );
    assert_eq!(
        fail(0, &stream(&[26, 1, 2, 3, 4, 5, 6, 7, 8])),
        Fail::Invalid
    );
    for token in (35u8..=42).chain(58..=60) {
        assert_eq!(
            fail(0, &stream(&[token])),
            Fail::Refuse(Rule::Xor),
            "{token}"
        );
    }
    // a super-block of length 0
    assert_eq!(fail(0, &stream(&[56, 0, 0, 0])), Fail::Invalid);
}

#[test]
fn sync_marks_check_the_block_since_the_bookmark() {
    // bookmark at block 0, 3 zero blocks, sync 3: matches; then an end
    let vector = ok(0, &stream(&[47, 9, 0, 3, 3, 50, 3, 19, 1, 0, 0]));
    assert_eq!(vector.ranges(), [((3 << 16) + 1, (3 << 16) + 1)]);
    assert_eq!(fail(0, &stream(&[47, 9, 0, 3, 3, 50, 2])), Fail::Invalid);
    // without a bookmark, the mark counts from block 0
    assert!(decode(0, &stream(&[3, 2, 51, 2, 0, 0])).is_ok());
    // 48-bit and 64-bit marks keep their low 32 bits
    assert!(decode(0, &stream(&[3, 1, 54, 1, 0, 0, 0, 7, 0, 0])).is_ok());
    assert!(decode(0, &stream(&[3, 1, 55, 1, 0, 0, 0, 7, 7, 7, 7, 0])).is_ok());
}

#[test]
fn the_interpolative_value_is_a_bijection_onto_its_range() {
    // `bic_value(r)` reads `logv` bits and maybe one more; every code gives a value of `0..=r`
    // and each value comes from one code.
    let ranges = (1u32..=1100).chain([65534, 65535, 65536, 1 << 20, (1 << 31) - 1, 1 << 31]);
    for r in ranges {
        let logv = 31 - (r + 1).leading_zeros();
        let (mut seen, mut count) = (std::collections::HashSet::new(), 0u64);
        let probe = |bits: u32| {
            let data = bits.to_le_bytes();
            let mut bytes = reader(&data);
            let mut bin = BitIn::new(&mut bytes);
            let x = bin.bic_value(r).unwrap();
            (x, bin.used)
        };
        if logv > 12 {
            // sample the codes
            for bits in (0u32..1 << 12).map(|k| k.wrapping_mul(2654435761) >> (31 - logv)) {
                let (x, _) = probe(bits);
                assert!(x <= r, "r {r} bits {bits:x} x {x}");
            }
            continue;
        }
        for bits in 0u32..1 << (logv + 1) {
            let (x, used) = probe(bits);
            assert!(x <= r, "r {r}");
            // each code once: a short code with the bit after it 0 only
            if used == logv && bits >> logv == 1 {
                continue;
            }
            assert!(seen.insert(x), "r {r} x {x} twice");
            count += 1;
        }
        assert_eq!(count, u64::from(r) + 1, "r {r}");
    }
    // `r + 1 == 0`: `bit_scan_reverse32(0)`
    let data = [0u8; 8];
    let mut bytes = reader(&data);
    assert_eq!(
        BitIn::new(&mut bytes).bic_value(u32::MAX),
        Err(Fail::Refuse(Rule::BicScan))
    );
}

#[test]
fn gamma_codes_and_their_undefined_lengths() {
    let mut out = BitOut::default();
    let values = [
        1u32,
        2,
        3,
        7,
        8,
        255,
        65535,
        65536,
        (1 << 31) - 1,
        1 << 31,
        u32::MAX,
    ];
    for &v in &values {
        out.gamma(v);
    }
    let data = out.bytes();
    let mut bytes = reader(&data);
    let mut bin = BitIn::new(&mut bytes);
    for &v in &values {
        assert_eq!(bin.gamma().unwrap(), v);
    }
    // 32 zero bits before the border bit
    let data = [0, 0, 0, 0, 1, 0, 0, 0];
    let mut bytes = reader(&data);
    assert_eq!(
        BitIn::new(&mut bytes).gamma(),
        Err(Fail::Refuse(Rule::GammaZeros))
    );
}

#[test]
fn a_gamma_gap_block() {
    // token 20: header (length 3, starts with 0), then the run ends 9 and 19 as gamma values
    // 10 (9 + 1) and 10
    let mut out = BitOut::default();
    out.gamma(10);
    out.gamma(10);
    let mut tokens = vec![20, 3 << 3, 0];
    tokens.extend(out.bytes());
    tokens.push(0);
    let vector = ok(0, &stream(&tokens));
    assert_eq!(vector.ranges(), [(10, 19)]);
    // length 0
    assert_eq!(fail(0, &stream(&[20, 0, 0])), Fail::Refuse(Rule::GapLength));
}

#[test]
fn gap_blocks_must_increase() {
    // token 14, length 3 starting with 1: runs end at 5, 9, 65535
    let vector = ok(0, &stream(&[14, 0x19, 0, 5, 0, 9, 0, 0]));
    assert_eq!(vector.ranges(), [(0, 5), (10, 65535)]);
    assert_eq!(vector.size(), 0);
    // a GAP block that is not one is kept; a query that reaches it runs the searches of NCBI's
    // three builds (scalar, SSE2, SSE4.2), answered where they agree (the answers of NCBI's three
    // builds of BitMagic 8.0.1), and counts as NCBI's `unsigned` sums do
    let vector = ok(0, &stream(&[14, 0x19, 0, 9, 0, 5, 0, 0]));
    for (n, bit, rank) in [
        (0, true, 1),
        (9, true, 10),
        (10, true, 15),
        (65535, true, 65540),
        (65536, false, 65540),
    ] {
        assert_eq!(vector.get_bit(n).unwrap(), bit, "{n}");
        assert_eq!(vector.count_to(n).unwrap(), rank, "{n}");
    }
    // an id list (token 18) that does not increase makes a GAP that does not either: bit 5 is 1
    // in the SSE4.2 build, 0 in the others
    let vector = ok(0, &stream(&[18, 2, 0, 9, 0, 5, 0, 0]));
    assert_eq!(vector.get_bit(5), Err(Fail::Refuse(Rule::GapSearch)));
    assert_eq!(vector.count_to(5).unwrap(), 0);
    assert!(vector.get_bit(9).unwrap());
    assert_eq!(vector.count_to(10).unwrap(), 2);
    // an empty id list reads `id_array_[0]`: memory that NCBI did not write, or the first id of
    // an earlier list (here block 0's, for block 1)
    assert_eq!(
        fail(0, &stream(&[18, 0, 0, 0])),
        Fail::Refuse(Rule::Unwritten)
    );
    let vector = ok(0, &stream(&[18, 1, 0, 5, 0, 18, 0, 0, 0]));
    assert_eq!(vector.ranges(), [(5, 5), (65541, 65541)]);
    let vector = ok(0, &stream(&[18, 3, 0, 5, 0, 6, 0, 9, 0, 0]));
    assert_eq!(vector.ranges(), [(5, 6), (9, 9)]);
    let ranks: Vec<u64> = [4, 5, 7, 8, 9, 10, 70000]
        .map(|n| vector.count_to(n).unwrap())
        .into();
    assert_eq!(ranks, [0, 1, 2, 2, 3, 3, 3]);
    // token 24 (ids off): an empty list is a full block
    let vector = ok(0, &stream(&[24, 0, 0, 0]));
    assert_eq!(vector.ranges(), [(0, 65535)]);
    let vector = ok(0, &stream(&[24, 3, 0, 0, 0, 1, 0, 0xFF, 0xFF, 0]));
    assert_eq!(vector.ranges(), [(2, 65534)]);
}

#[test]
fn only_a_query_that_reaches_a_block_that_is_not_one_is_refused() {
    // block 3 is a GAP block that is not one (65535 is not its last entry): the rows of blocks 0
    // to 2 are answered, as NCBI does not read it for them; in block 3 its answer splits NCBI's
    // builds (bits 196708.. are 1 in the SSE4.2 build only) where it is refused, and the ranks
    // past it are NCBI's `unsigned` sums
    let data = hex("020180000001000200050000080013050007020000000E2000409C0A0050C300");
    let vector = ok(524288, &data);
    for (n, bit, rank) in [
        (0, false, 0),
        (5, true, 1),
        (6, false, 1),
        (65536, false, 1),
        (196607, false, 1),
        (196608, false, 1),
        (262143, true, 4294942842),
        (262144, false, 4294942842),
    ] {
        assert_eq!(vector.get_bit(n).unwrap(), bit, "{n}");
        assert_eq!(vector.count_to(n).unwrap(), rank, "{n}");
    }
    for n in [196708, 196709, 196808] {
        assert_eq!(vector.get_bit(n), Err(Fail::Refuse(Rule::GapSearch)), "{n}");
        assert_eq!(vector.count_to(n).unwrap(), 1, "{n}");
    }
    // token 61's exceptions set bits in a GAP block that is not one (ends 60100, 60201, 54765,
    // 54867): `gap_set_value` of NCBI's three builds leaves it the same, so it is decoded, and
    // its rows are NCBI's
    let data = hex("02018000000100020005000000803DB050C800DEFC030082A915603C9E39000000000000");
    let vector = ok(2147483648, &data);
    assert!(!vector.get_bit(5).unwrap());
    assert_eq!(vector.count_to(5).unwrap(), 0);
}

#[test]
fn bits_set_in_a_gap_block_that_is_not_one() {
    // t2e: a GAP block that is not one (token 20), then a bit set in it (token 19): NCBI's three
    // builds' `gap_set_value` leave the block the same, and its rows are NCBI's; bit 100 is 1 in
    // the SSE4.2 build only, and the ranks are NCBI's `unsigned` sums
    let data = hex("0201800000010002000500000080141800C01200F0FF0F000007FFFFFFFF13050000");
    let vector = ok(2147483648, &data);
    for (n, bit, rank) in [
        (0, false, 0),
        (4, false, 0),
        (5, true, 1),
        (9, false, 1),
        (99, false, 1),
        (4660, false, 0),
        (65536, false, 0),
    ] {
        assert_eq!(vector.get_bit(n).unwrap(), bit, "{n}");
        assert_eq!(vector.count_to(n).unwrap(), rank, "{n}");
    }
    assert_eq!(vector.get_bit(100), Err(Fail::Refuse(Rule::GapSearch)));
    assert_eq!(vector.count_to(100).unwrap(), 1);
    // a bit set where the builds differ: the block of the id list 9, 5, whose bit 5 is 0 in the
    // scalar and SSE2 builds (which set it) and 1 in the SSE4.2 build (which leaves the block)
    let vector = ok(0, &stream(&[18, 2, 0, 9, 0, 5, 0, 0]));
    let segs = &vector.one().expect("one vector").segs;
    let Kind::Gap(entries) = &segs[0].kind else {
        panic!("{:?}", segs[0].kind);
    };
    let gap = Arr::of_block(
        entries.to_vec(),
        GAP_LEVELS[usize::from(gap_level(entries[0]))],
    );
    assert_eq!(
        gap_set_value(&mut gap.clone(), true, 5, false, None),
        Err(Fail::Refuse(Rule::GapSearch))
    );
    // each build alone: the scalar and SSE2 builds set the bit (a new run 5..=5), the SSE4.2
    // build finds it set and leaves the block
    for (build, changed) in [
        (Build::Scalar, true),
        (Build::Sse2, true),
        (Build::Sse42, false),
    ] {
        let mut one = gap.clone();
        gap_set_value(&mut one, true, 5, false, Some(build)).unwrap();
        assert_eq!(one.written() != gap.written(), changed, "{build:?}");
    }
}

#[test]
fn a_count_past_broken_blocks_adds_each_builds_sums() {
    // review 2's c02399: blocks whose scalar and SSE whole-block counts differ by multiples of
    // 65536 in opposite ways, so that past them each build's `unsigned` total is the same
    let data = hex(concat!(
        "180EA10045156B1A301DD1336A41B641EF49185D778498896D92D89A57AA6DB71CDC66078DEDE8ED5EEF02",
        "0E9800FAC2FAF37D1BC17C9A426F920100600A1F4E78B0677EB2CA8EF3830E6FAFF8CB45E537F10200"
    ));
    let vector = ok(u32::MAX, &data);
    assert_eq!(vector.count_to(196608).unwrap(), 219241);
}

#[test]
fn bits_that_leave_the_builds_different_blocks() {
    // review 2's E1, E2, E0: token 18 makes block 0 the GAP block of the id list 9, 5 (not one);
    // token 19 sets bit 5 in it, which the scalar and SSE2 builds' `gap_set_value` add and the
    // SSE4.2 build's finds set: the data is decoded once for each build
    let rows = |vector: &BitVector, rows: &[u32]| -> Vec<String> {
        let answer = |result: Result<u64, Fail>| match result {
            Ok(value) => value.to_string(),
            Err(_) => "?".into(),
        };
        rows.iter()
            .map(|&n| {
                let bit = answer(vector.get_bit(n).map(u64::from));
                format!("{n}:{bit}:{}", answer(vector.count_to(n)))
            })
            .collect()
    };
    let rows_e = [0, 4, 5, 6, 9, 10, 65535, 65536, 65540, 131071];
    // E1: block 0 then made full (token 2): all builds hold the same vector
    let e1 = ok(
        u32::MAX,
        &hex("181202000900050007FFFFFFFF13050007FFFFFFFF0200"),
    );
    assert!(e1.one().is_some());
    assert_eq!(
        rows(&e1, &rows_e).join(","),
        "0:1:1,4:1:5,5:1:6,6:1:7,9:1:10,10:1:11,65535:1:65536,65536:0:65536,65540:0:65536,\
         131071:0:65536"
    );
    // E2: block 1 gets bits too; the bits are the same in all builds, the ranks from row 5 on
    // differ (scalar and SSE2 65539 at 131071, SSE4.2 65538)
    let e2 = ok(u32::MAX, &hex("181202000900050007FFFFFFFF1305000200"));
    assert!(e2.one().is_none());
    assert_eq!(
        rows(&e2, &rows_e).join(","),
        "0:0:0,4:0:0,5:1:?,6:0:?,9:1:?,10:0:?,65535:0:?,65536:1:?,65540:1:?,131071:1:?"
    );
    // E0: the bit is not set; bit 5 is 1 in the SSE4.2 build only
    let e0 = ok(u32::MAX, &hex("181202000900050000"));
    assert_eq!(
        rows(&e0, &rows_e[..7]).join(","),
        "0:0:0,4:0:0,5:?:0,6:0:0,9:1:1,10:0:2,65535:0:2"
    );
}

#[test]
fn bits_set_in_a_gap_block_that_is_one_need_one_search() {
    // the bits of a GAP block that is one, set one by one in any order (token 19, then token 7
    // back to block 0), give the block of those bits, also past a level's limit
    // (`extend_gap_block`) and past 1276 entries (a bitmap)
    let mut state = 12345u32;
    let mut next = || {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        state >> 8
    };
    for count in [1usize, 2, 3, 60, 61, 62, 63, 130, 700, 1000] {
        let mut tokens = Vec::new();
        let mut bits = std::collections::BTreeSet::new();
        for _ in 0..count {
            let bit = (next() % 3000 * 7) as u16;
            bits.insert(bit);
            tokens.push(19);
            tokens.extend(bit.to_le_bytes());
            tokens.push(7);
            tokens.extend(u32::MAX.to_le_bytes());
        }
        tokens.push(0);
        let vector = ok(0, &stream(&tokens));
        let want: Vec<(u32, u32)> = bits.iter().map(|&b| (u32::from(b), u32::from(b))).collect();
        let got = vector.ranges();
        let merged: Vec<(u32, u32)> = want.iter().fold(Vec::new(), |mut out, &(a, b)| {
            match out.last_mut() {
                Some(last) if last.1 + 1 == a => last.1 = b,
                _ => out.push((a, b)),
            }
            out
        });
        assert_eq!(got, merged, "{count}");
    }
}

#[test]
fn long_gaps_are_converted_to_bits() {
    // token 14 with 1300 entries (above 1276): converted, so only its runs of 1s count
    let mut ends = vec![10u16, 20, 5, 7];
    ends.extend((5..1300).map(|k| 100 + k as u16));
    let mut tokens = vec![14];
    tokens.extend((1300u16 << 3).to_le_bytes());
    for e in &ends {
        tokens.extend(e.to_le_bytes());
    }
    tokens.push(0);
    let vector = ok(0, &stream(&tokens));
    for (n, bit) in [
        (5, false),
        (6, true),
        (7, true),
        (10, false),
        (11, true),
        (20, true),
    ] {
        assert_eq!(vector.get_bit(n).unwrap(), bit, "{n}");
    }
    for (n, bit) in [
        (21, false),
        (106, true),
        (107, false),
        (1398, true),
        (1399, false),
    ] {
        assert_eq!(vector.get_bit(n).unwrap(), bit, "{n}");
    }
    assert_eq!(vector.count(), 2 + 10 + 647 + (65535 - 1400 + 1));
    // a run of 1s from 41 back to 32 or 31: `or_bit_block(dest, 41, 2^32 - d)` writes nothing
    // while `d` is below 41's bit in its word (9), shifts by 32 at 9 and writes past the block
    // above
    let gap14 = |head: &[u16]| {
        let mut tokens = vec![14];
        tokens.extend((1300u16 << 3).to_le_bytes());
        for k in 0..1299 {
            let end = head.get(k).copied().unwrap_or(100 + k as u16);
            tokens.extend(end.to_le_bytes());
        }
        tokens.push(0);
        stream(&tokens)
    };
    let vector = ok(0, &gap14(&[40, 32, 50, 60]));
    for (n, bit) in [
        (31, false),
        (41, false),
        (50, false),
        (51, true),
        (60, true),
        (61, false),
    ] {
        assert_eq!(vector.get_bit(n).unwrap(), bit, "{n}");
    }
    assert_eq!(
        fail(0, &gap14(&[40, 31, 50, 60])),
        Fail::Refuse(Rule::GapOrder)
    );
    // with 1275 entries it stays a GAP block, which must increase
    let mut tokens = vec![14];
    tokens.extend((1275u16 << 3).to_le_bytes());
    for e in &ends[..1274] {
        tokens.extend(e.to_le_bytes());
    }
    tokens.push(0);
    // with 1275 entries it stays a GAP block that is not one: bits 6 and 7 are 1 in the scalar and
    // SSE4.2 builds, 11 and 20 in the SSE2 and SSE4.2 builds
    let vector = ok(0, &stream(&tokens));
    for n in [6, 7, 11, 20] {
        assert_eq!(vector.get_bit(n), Err(Fail::Refuse(Rule::GapSearch)), "{n}");
    }
    for (n, bit, rank) in [
        (5, false, 0),
        (21, false, 12),
        (106, true, 13),
        (200, true, 60),
        (65535, false, 647),
    ] {
        assert_eq!(vector.get_bit(n).unwrap(), bit, "{n}");
        assert_eq!(vector.count_to(n).unwrap(), rank, "{n}");
    }
    assert_eq!(vector.count_to(20).unwrap(), 10);
    // token 18 with 4097 ids 1, 3, ..., 8193: a GAP of 8195 entries, whose header keeps 3; the
    // conversion reads the runs to entry 3 only
    let mut tokens = vec![18];
    tokens.extend(4097u16.to_le_bytes());
    for k in 0..4097u16 {
        tokens.extend((2 * k + 1).to_le_bytes());
    }
    tokens.push(0);
    let vector = ok(0, &stream(&tokens));
    assert_eq!(vector.ranges(), [(1, 1)]);
}

#[test]
fn interpolative_tokens() {
    // token 31: min 3, max 900, 5 values with 100, 200, 300 between
    let mut out = BitOut::default();
    out.bic(&[100, 200, 300], 3, 900);
    let mut tokens = vec![31, 3, 0, 0x84, 3, 5, 0];
    tokens.extend(out.bytes());
    tokens.push(0);
    let vector = ok(0, &stream(&tokens));
    assert_eq!(
        vector.ranges(),
        [(3, 3), (100, 100), (200, 200), (300, 300), (900, 900)]
    );
    // token 28 (GAP id list): min 10, max 20, gamma(len - 4) with len 5
    let mut out = BitOut::default();
    out.gamma(1);
    out.bic(&[11, 12, 15], 10, 20);
    let mut tokens = vec![28, 10, 0, 20, 0];
    tokens.extend(out.bytes());
    tokens.push(0);
    let vector = ok(0, &stream(&tokens));
    assert_eq!(vector.ranges(), [(10, 12), (15, 15), (20, 20)]);
    // a super-block (56): index 1, 3 values from min 5 to 16777215 - 2
    let mut out = BitOut::default();
    out.bic(&[70000], 5, (1 << 24) - 2);
    let mut tokens = vec![56, 0x00, 1, 3, 5, 2];
    tokens.extend(out.bytes());
    tokens.push(0);
    let vector = ok(0, &stream(&tokens));
    let base = 1u32 << 24;
    assert_eq!(
        vector.ranges(),
        [
            (base + 5, base + 5),
            (base + 70000, base + 70000),
            (base + (1 << 24) - 2, base + (1 << 24) - 2)
        ]
    );
    // the size stays
    assert_eq!(vector.size(), 0);
}

#[test]
fn bit_blocks() {
    // token 11: 2048 words, word 1 = 0x80000001
    let mut tokens = vec![3, 2, 11];
    let mut words = vec![0u8; 8192];
    words[4..8].copy_from_slice(&0x8000_0001u32.to_le_bytes());
    tokens.extend(&words);
    tokens.push(0);
    let vector = ok(0, &stream(&tokens));
    let at = 2 << 16;
    assert_eq!(vector.ranges(), [(at + 32, at + 32), (at + 63, at + 63)]);
    // token 17: words 2..=3
    let vector = ok(0, &stream(&[17, 2, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0x80, 0]));
    assert_eq!(vector.ranges(), [(64, 64), (127, 127)]);
    // `tail == head - 1`: no word (little-endian), the block made empty
    let vector = ok(0, &stream(&[17, 3, 0, 2, 0, 0]));
    assert_eq!(vector.ranges(), []);
    assert_eq!(
        fail(0, &stream(&[17, 5, 0, 2, 0, 0])),
        Fail::Refuse(Rule::Interval)
    );
    // into a null block the first loop clears `head` words, past the block above 2048
    assert_eq!(
        fail(0, &stream(&[17, 1, 8, 0, 8, 0])),
        Fail::Refuse(Rule::Interval)
    );
    // token 22: a run of 1 zero word, 1 word, 2046 zero words
    let vector = ok(
        0,
        &stream(&[22, 0, 1, 0, 1, 0, 0xFF, 0, 0, 0, 0xFE, 0x07, 0]),
    );
    assert_eq!(vector.ranges(), [(32, 39)]);
    assert_eq!(
        fail(0, &stream(&[22, 1, 0x01, 0x08])),
        Fail::Refuse(Rule::ZeroRuns)
    );
    // token 30: all bits but 0 and 65535
    let vector = ok(0, &stream(&[30, 2, 0, 0, 0, 0xFF, 0xFF, 0]));
    assert_eq!(vector.ranges(), [(1, 65534)]);
    // token 34: digest wave 1 (words 32..64), first word 1
    let mut tokens = vec![34, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];
    tokens.extend([0u8; 124]);
    tokens.push(0);
    let vector = ok(0, &stream(&tokens));
    assert_eq!(vector.ranges(), [(1024, 1024)]);
}

#[test]
fn queries_match_a_plain_bitmap() {
    // blocks 0, 1 and 3 from words of a fixed generator, block 2 full
    let mut state = 0x1234_5678_9abc_def0u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut tokens = vec![];
    let mut bits: Vec<bool> = vec![];
    for nb in 0..4 {
        if nb == 2 {
            tokens.push(2);
            bits.extend(std::iter::repeat_n(true, 65536));
            continue;
        }
        tokens.push(11);
        for _ in 0..2048 {
            let r = next();
            let word = match r % 4 {
                0 => 0,
                1 => u32::MAX,
                _ => (r >> 32) as u32,
            };
            tokens.extend(word.to_le_bytes());
            bits.extend((0..32).map(|k| word >> k & 1 != 0));
        }
    }
    tokens.push(0);
    let vector = ok(4 << 16, &stream(&tokens));
    let mut rank = 0u64;
    for (n, &bit) in bits.iter().enumerate() {
        rank += u64::from(bit);
        if n % 97 == 0 || n % 65536 < 3 || n % 65536 > 65533 {
            assert_eq!(vector.get_bit(n as u32).unwrap(), bit, "{n}");
            assert_eq!(vector.count_to(n as u32).unwrap(), rank, "{n}");
        }
    }
    assert_eq!(vector.count(), rank);
    assert_eq!(vector.count_to(u32::MAX).unwrap(), rank);
    assert!(!vector.get_bit(4 << 16).unwrap());
}

#[test]
fn big_endian_streams() {
    // bo 0: 16-bit values big-endian; the 16-bit array repeats its first value
    let vector = ok(0, &[0x10, 0x00, 19, 0x01, 0x02, 0]);
    assert_eq!(vector.ranges(), [(0x102, 0x102)]);
    // token 16 (bits on) with 2 ids: read one by one
    let vector = ok(0, &[0x10, 0x00, 16, 0, 2, 0, 5, 0, 9, 0]);
    assert_eq!(vector.ranges(), [(5, 5), (9, 9)]);
    // token 18 with 2 ids: both are the first one, a GAP that does not increase, whose bit 5 is 0
    // in the SSE4.2 build and 1 in the others
    let vector = ok(0, &[0x10, 0x00, 18, 0, 2, 0, 5, 0, 9, 0]);
    assert_eq!(vector.get_bit(5), Err(Fail::Refuse(Rule::GapSearch)));
    assert!(!vector.get_bit(6).unwrap());
    assert_eq!(vector.count_to(5).unwrap(), 1);
    assert_eq!(vector.count_to(6).unwrap(), 2);
    // ... with 1 id: one value, two bytes
    let vector = ok(0, &[0x10, 0x00, 18, 0, 1, 0, 5, 0]);
    assert_eq!(vector.ranges(), [(5, 5)]);
    // a GAP block of length 1 into a null block reads nothing; into a block that exists the
    // big-endian `get_16(p, 0)` reads two bytes (the little-endian one none), and the GAP of all
    // ones ORed into a GAP block makes it the full block
    let vector = ok(0, &[0x10, 0x00, 14, 0x00, 0x09, 0]);
    assert_eq!(vector.ranges(), [(0, 65535)]);
    let vector = ok(
        0,
        &[
            0x10, 0x00, 19, 0, 1, 7, 0xFF, 0xFF, 0xFF, 0xFF, 14, 0x00, 0x09, 0xAB, 0xCD, 0,
        ],
    );
    assert_eq!(vector.ranges(), [(0, 65535)]);
    let vector = ok(
        0,
        &stream(&[19, 1, 0, 7, 0xFF, 0xFF, 0xFF, 0xFF, 14, 0x09, 0x00, 0]),
    );
    assert_eq!(vector.ranges(), [(0, 65535)]);
}

/// The differential driver: each line of `BM_DIFF_IN` is `ID SIZE HEX ROWS` (ROWS a comma list
/// or `-`); each line of `BM_DIFF_OUT` is `ID ok SIZE COUNT CONSUMED HASH RUNS QUERIES`,
/// `ID invalid` or `ID refuse RULE`, where HASH is FNV-1a 64 over the runs of set bits (first
/// and last as little-endian 32-bit values) and each query is `ROW:BIT:RANK` (`get_bit`,
/// `count_range(0, ROW)` as NCBI's 32-bit `size_type`).
#[test]
fn differential() {
    let Ok(input) = std::env::var("BM_DIFF_IN") else {
        return;
    };
    let output = std::env::var("BM_DIFF_OUT").expect("BM_DIFF_OUT");
    let text = std::fs::read_to_string(input).unwrap();
    let mut out = String::new();
    for line in text.lines() {
        let mut fields = line.split(' ');
        let (Some(id), Some(size), Some(data), Some(rows)) =
            (fields.next(), fields.next(), fields.next(), fields.next())
        else {
            continue;
        };
        let size = size.parse::<u32>().unwrap();
        let data = hex(data);
        match decode_to_end(size, &data) {
            Err(Fail::Invalid) => out.push_str(&format!("{id} invalid\n")),
            Err(Fail::Refuse(rule)) => out.push_str(&format!("{id} refuse {rule:?}\n")),
            Ok((vector, end)) => {
                // the harness's walk of NCBI's blocks; `-` where a block's entries are unknown
                let (count, hash, runs) = match harness_walk(&vector) {
                    Some(runs) => {
                        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
                        for &(a, b) in &runs {
                            for byte in a.to_le_bytes().into_iter().chain(b.to_le_bytes()) {
                                hash ^= u64::from(byte);
                                hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
                            }
                        }
                        let count: u64 = runs.iter().map(|&(a, b)| u64::from(b - a) + 1).sum();
                        (
                            count.to_string(),
                            format!("{hash:016x}"),
                            runs.len().to_string(),
                        )
                    }
                    None => ("-".into(), "-".into(), "-".into()),
                };
                let answer = |result: Result<u64, Fail>| match result {
                    Ok(value) => value.to_string(),
                    Err(Fail::Refuse(rule)) => format!("?{rule:?}"),
                    Err(Fail::Invalid) => "?Invalid".into(),
                };
                let queries: Vec<String> = rows
                    .split(',')
                    .filter(|row| *row != "-")
                    .map(|row| {
                        let row = row.parse::<u32>().unwrap();
                        // `count_range` counts to `bm::id_max - 1` at most (bm.h:3251-3252)
                        format!(
                            "{row}:{}:{}",
                            answer(vector.get_bit(row).map(u64::from)),
                            answer(
                                vector
                                    .count_to(row.min(u32::MAX - 1))
                                    .map(|n| n as u32 as u64)
                            )
                        )
                    })
                    .collect();
                let queries = if queries.is_empty() {
                    "-".to_string()
                } else {
                    queries.join(",")
                };
                out.push_str(&format!(
                    "{id} ok {} {count} {end} {hash} {runs} {queries}\n",
                    vector.size(),
                ));
            }
        }
    }
    std::fs::write(output, out).unwrap();
}

/// The runs that the differential harness's walk of NCBI's blocks reports (`walk`, merging a run
/// that starts right after the last): a GAP block's runs of 1s whose end is not before their start,
/// up to its length or 65535; `None` for a block whose entries the data does not decide.
fn harness_walk(vector: &BitVector) -> Option<Vec<(u32, u32)>> {
    let mut out: Vec<(u32, u32)> = Vec::new();
    let mut put = |a: u32, b: u32| match out.last_mut() {
        Some(last) if last.1.checked_add(1) == Some(a) => last.1 = b,
        _ => out.push((a, b)),
    };
    for seg in &vector.one()?.segs {
        match &seg.kind {
            Kind::Poison(_) => return None,
            Kind::Gap(gap) => {
                let base = seg.first << 16;
                let mut value = gap[0] & 1 != 0;
                let mut from = 0u32;
                let ends = gap.get(1..=usize::from(gap[0] >> 3))?;
                for &end in ends {
                    let end = u32::from(end);
                    if value && end >= from {
                        put(base + from, base + end);
                    }
                    from = end + 1;
                    value = !value;
                    if end >= 65535 {
                        break;
                    }
                }
            }
            _ => seg.each_run(&mut put),
        }
    }
    Some(out)
}
