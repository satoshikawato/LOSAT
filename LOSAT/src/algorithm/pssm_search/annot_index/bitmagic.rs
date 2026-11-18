//! BitMagic 8.0.1's deserialization as NCBI BLAST+ 2.17.0's object manager runs it on the data of
//! a Seq-table's `BVector-data` (line RG Q-BM, `bm/SPEC.md`): `bm::deserialize` into a
//! `bm::bvector<>` made with the `size` member, in NCBI's 32-bit BitMagic build (`BM64ADDR` not
//! defined; `BM_ASSERT` compiled out). NCBI passes the data without its length; LOSAT decodes the
//! same tokens in NCBI's order into the blocks that NCBI's vector holds (none, a GAP block, a
//! bitmap, or the full block's pointer), by NCBI's steps for each, and stops at the first event of
//! NCBI's that it does not port (`Fail::Refuse`): a read past the data, an XOR reference token, a
//! step that NCBI runs with undefined behaviour or a write past its buffers, a read of memory that
//! NCBI did not write in this decoding, or a search of a GAP block that is not one whose result
//! differs between NCBI's builds. NCBI's own checks throw
//! `std::logic_error("BM::Invalid serialization format")` (`Fail::Invalid`).
//!
//! The decoded vector keeps the size and the set bits, which is all that NCBI's index asks of it
//! (`BitVector`: the size, a bit, and the count of set bits up to a bit; `bm/SPEC.md` §1.3). A GAP
//! block that is not one (its ends do not increase up to 65535) is kept as NCBI keeps it. A step
//! that searches it, a bit that the decoding sets in it (`gap_set_value`) or a query that reaches
//! it, runs the searches of NCBI's three builds of BitMagic (scalar, SSE2, SSE4.2: which one NCBI
//! runs depends on the machine it was built for) and goes on where the three agree, refused
//! (`Rule::GapSearch`) where they differ. Rows that do not reach it are answered.
//!
//! NCBI reference (598d8ae6): c++/src/objects/seqtable/BVector_data.cpp:73-82
//! ```c++
//! void CBVector_data::x_CreateBitVector(void) const
//! {
//!     CMutexGuard guard(sx_PrepareMutex_bvector);
//!     if ( !m_BitVector ) {
//!         // make bvector
//!         AutoPtr<bm::bvector<> > bv(new bm::bvector<>(GetSize()));
//!         bm::deserialize(*bv, (const unsigned char*)&GetData()[0]);
//!         m_BitVector.reset(bv.release());
//!     }
//! }
//! ```
//! The module has no `crate::` imports, so that `rustc --edition 2021 --test bitmagic.rs` builds it
//! and its tests alone (the differential driver of `bitmagic_tests.rs`).

/// Blocks of a vector (`bm::set_total_blocks`) and bits of a block (`bm::gap_max_bits`).
const TOTAL_BLOCKS: u32 = 65536;
const BLOCK_BITS: u32 = 65536;
/// A block of 65536 bits as 64-bit words (bit `p` is bit `p % 64` of word `p / 64`, the same
/// numbering as BitMagic's 32-bit words).
type Block = [u64; 1024];
/// The entries of `gap_temp_block_`, `id_array_` and `sb_id_array_` (bmserial.h:5122-5128) and
/// of `ex0_arr_` / `ex1_arr_` (`alloc_bit_block(tmp_buff_alloc_factor)`, 8 blocks of 8192 bytes
/// as `gap_word_t`, bmserial.h:1269, 5119-5120).
const GAP_ENTRIES: usize = 65536;
const EX_ENTRIES: usize = 32768;
/// The entries of the GAP result buffer of `combine_operation_block_or` and
/// `combine_operation_with_block` (`tmp_buf[bm::gap_equiv_len * 3]`, bm.h:7287, 7665).
const OR_BUF: usize = 3 * 4096;
/// The decoded vector keeps a block of more runs than this as its words.
const RUNS_MAX: usize = 2048;

/// Why the data is not decoded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Fail {
    /// NCBI's `std::logic_error("BM::Invalid serialization format")` (bmserial.h:620): the header
    /// flag of a 64-bit vector, the tokens 25 and 26, a sync mark that does not match, a
    /// super-block of length 0, an unknown token.
    Invalid,
    /// An event that LOSAT does not port (`Rule`).
    Refuse(Rule),
}

/// The event that refused the data (`bm/SPEC.md` §3.2, §3.3; LOSAT's tests and differential
/// driver read it).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Rule {
    /// A byte or a bit-stream word past the data (U1).
    ReadPastEnd,
    /// An XOR reference token, which dereferences NCBI's null reference vector (U2).
    Xor,
    /// 32 or more zero bits before an Elias gamma value (U3).
    GammaZeros,
    /// `bit_scan_reverse32(0)` in the interpolative decoder (U4).
    BicScan,
    /// A value of the interpolative decoder past the 65536 bits of a block (U4).
    BicWrite,
    /// An index past one of NCBI's arrays (U4, U7, U8).
    Bounds,
    /// A bit-interval block whose words are not within the block (U5).
    Interval,
    /// A zero-runs block whose words run past the block (U6).
    ZeroRuns,
    /// A GAP block of length 0 that NCBI reads (U7, U9), or one whose length left the header's
    /// 13 bits before NCBI searches it.
    GapLength,
    /// A run of a GAP block that is not one (its ends do not increase up to 65535) that ends
    /// before it starts where NCBI's `or_bit_block` shifts by 32 or writes past the block (U12).
    GapOrder,
    /// An entry of NCBI's work arrays or of a GAP block that this decoding did not write.
    Unwritten,
    /// A window of 0 words: a division by zero (U11).
    Window,
    /// A run of bits past the block (U10).
    RunRange,
    /// A search of a GAP block that is not one (a query, or `gap_set_value` while decoding) whose
    /// result differs between NCBI's builds (scalar, SSE2, SSE4.2): NCBI's result depends on the
    /// machine it was built for.
    GapSearch,
}

/// The `what()` of NCBI's `std::logic_error` (`Fail::Invalid`).
pub(crate) const INVALID: &str = "BM::Invalid serialization format";

const fn refuse(rule: Rule) -> Fail {
    Fail::Refuse(rule)
}

/// The data's size and set bits as NCBI's index asks them (`bm::bvector<>::size`, `get_bit`,
/// `count_range`), kept by blocks: runs of full blocks, blocks as runs of set bits, blocks as
/// bitmaps, and the blocks that NCBI keeps as GAP blocks that are not one, which a query that
/// reaches them searches as NCBI's builds do.
#[derive(Debug)]
pub(crate) struct BitVector {
    size: u32,
    segs: Vec<Seg>,
    count: u64,
    /// The whole-block counts of the blocks that are not one, added as NCBI adds them
    /// (`unsigned`), or the refusal of the first that has none.
    other: Result<u32, Fail>,
    /// Whether the vector has such a block (its counts are then NCBI's `unsigned` sums).
    any_other: bool,
}

#[derive(Debug)]
struct Seg {
    first: u32,
    last: u32,
    /// The set bits of the segments before, but the blocks that are not one.
    before: u64,
    /// `BitVector::other` of the segments before.
    other: Result<u32, Fail>,
    kind: Kind,
}

#[derive(Debug)]
enum Kind {
    Full,
    /// Runs `(first, last)` in order, and the set bits of the runs before each (one more entry:
    /// all of them).
    Runs(Box<[(u16, u16)]>, Box<[u32]>),
    /// The words, and the set bits of the 16-word chunks before each.
    Bits(Box<Block>, Box<[u32; 64]>),
    /// A GAP block that is not one: the entries of its buffer that NCBI wrote.
    Gap(Box<[u16]>),
    /// A block whose entries depend on memory or code that the data does not decide: refused by
    /// any query that reaches it.
    Poison(Rule),
}

impl BitVector {
    /// `bvector<>::size()`.
    pub(crate) fn size(&self) -> u32 {
        self.size
    }

    fn seg(&self, block: u32) -> (usize, Option<&Seg>) {
        let at = self.segs.partition_point(|seg| seg.last < block);
        (at, self.segs.get(at).filter(|seg| seg.first <= block))
    }

    /// `bvector<>::get_bit(n)` (bm.h:3604-3631): a GAP block's `gap_test_unr`.
    pub(crate) fn get_bit(&self, n: u32) -> Result<bool, Fail> {
        let p = n & (BLOCK_BITS - 1);
        Ok(match self.seg(n >> 16).1.map(|seg| &seg.kind) {
            None => false,
            Some(Kind::Full) => true,
            Some(Kind::Runs(runs, _)) => {
                let at = runs.partition_point(|run| u32::from(run.1) < p);
                runs.get(at).is_some_and(|run| u32::from(run.0) <= p)
            }
            Some(Kind::Bits(words, _)) => words[(p >> 6) as usize & 1023] >> (p & 63) & 1 != 0,
            Some(Kind::Gap(gap)) => gap_test3(gap, p)?,
            Some(Kind::Poison(rule)) => return Err(refuse(*rule)),
        })
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bm.h:3258-3330
    // (`count_range_no_check(0, right)`: the first block whole by `block_bitcount` when `right`
    // is past it, else by `gap_bit_count_range`; the blocks between by `block_bitcount`; the
    // last by `gap_bit_count_to`; summed as `unsigned`)
    /// `bvector<>::count_range(0, n)` (bm.h:3244-3254): the set bits of `[0, n]`.
    pub(crate) fn count_to(&self, n: u32) -> Result<u64, Fail> {
        let block = n >> 16;
        let p = n & (BLOCK_BITS - 1);
        let (at, seg) = self.seg(block);
        let total = |valid: u64, other: Result<u32, Fail>, inside: u32| -> Result<u64, Fail> {
            if !self.any_other {
                return Ok(valid);
            }
            Ok(u64::from(
                (valid as u32).wrapping_add(other?).wrapping_add(inside),
            ))
        };
        let Some(seg) = seg else {
            return match self.segs.get(at) {
                Some(seg) => total(seg.before, seg.other, 0),
                None => total(self.count, self.other, 0),
            };
        };
        let inside = match &seg.kind {
            Kind::Full => u64::from(block - seg.first) * u64::from(BLOCK_BITS) + u64::from(p) + 1,
            Kind::Runs(runs, before) => {
                let i = runs.partition_point(|run| u32::from(run.1) < p);
                let partial = runs
                    .get(i)
                    .filter(|run| u32::from(run.0) <= p)
                    .map_or(0, |run| p - u32::from(run.0) + 1);
                u64::from(before[i] + partial)
            }
            Kind::Bits(words, chunks) => {
                let w = (p >> 6) as usize;
                let chunk = w >> 4;
                let mut sum = chunks[chunk & 63];
                for word in words.iter().take(w).skip(chunk << 4) {
                    sum += word.count_ones();
                }
                let mask = u64::MAX >> (63 - (p & 63));
                sum += (words[w & 1023] & mask).count_ones();
                u64::from(sum)
            }
            Kind::Gap(gap) => {
                let inside = if block != 0 {
                    gap_count_to(gap, p)?
                } else if p == BLOCK_BITS - 1 {
                    gap_count3(gap)?
                } else {
                    gap_count_range(gap, 0, p)?
                };
                return total(seg.before, seg.other, inside);
            }
            Kind::Poison(rule) => return Err(refuse(*rule)),
        };
        total(seg.before + inside, seg.other, 0)
    }

    /// All set bits.
    #[cfg(test)]
    pub(crate) fn count(&self) -> u64 {
        self.count
    }

    /// The runs of set bits `(first, last)` in order, of a vector without blocks that are not
    /// one.
    #[cfg(test)]
    pub(crate) fn ranges(&self) -> Vec<(u32, u32)> {
        let mut out: Vec<(u32, u32)> = Vec::new();
        let mut push = |a: u32, b: u32| match out.last_mut() {
            Some(last) if last.1.checked_add(1) == Some(a) => last.1 = b,
            _ => out.push((a, b)),
        };
        for seg in &self.segs {
            for (a, b) in seg.runs() {
                push(a, b);
            }
        }
        out
    }
}

impl Seg {
    /// The runs of set bits of a segment, in order (none for a block that is not one).
    #[cfg(test)]
    fn runs(&self) -> Vec<(u32, u32)> {
        let base = self.first << 16;
        match &self.kind {
            Kind::Full => vec![(base, (self.last << 16) | 0xFFFF)],
            Kind::Runs(runs, _) => runs
                .iter()
                .map(|&(a, b)| (base | u32::from(a), base | u32::from(b)))
                .collect(),
            Kind::Bits(words, _) => block_runs(words, usize::MAX)
                .unwrap_or_default()
                .into_iter()
                .map(|(a, b)| (base | u32::from(a), base | u32::from(b)))
                .collect(),
            Kind::Gap(_) | Kind::Poison(_) => Vec::new(),
        }
    }
}

/// `bm::deserialize(bv, data)` (bmserial.h:4085-4125) of a `bm::bvector<>(size)`
/// (BVector_data.cpp:77-78).
pub(crate) fn decode(size: u32, data: &[u8]) -> Result<BitVector, Fail> {
    decode_to_end(size, data).map(|(vector, _)| vector)
}

/// `decode`, and the value `bm::deserialize` returns (`dec.size()`, the bytes it moved over).
pub(crate) fn decode_to_end(size: u32, data: &[u8]) -> Result<(BitVector, usize), Fail> {
    let mut decoder = Decoder {
        bytes: Bytes {
            data,
            pos: 0,
            beq: false,
        },
        size,
        store: Store::default(),
        s: Scratch::default(),
    };
    let end = decoder.run()?;
    Ok((decoder.store.finish(decoder.size)?, end))
}

// ------------------------------------------------------------------------------------------------
// The byte decoders.

/// `bm::decoder` (little-endian, encoding.h:919-1152) or, for a byte-order byte other than 1,
/// `bm::decoder_little_endian` (encoding.h:1158-1295: big-endian 16- and 32-bit values, and
/// little-endian 24- and 48-bit ones, NCBI's own 64-bit order, the 16-bit array that repeats its
/// first value). `pos` may pass the data's end without a read (that array).
struct Bytes<'a> {
    data: &'a [u8],
    pos: usize,
    beq: bool,
}

impl Bytes<'_> {
    fn take<const N: usize>(&mut self) -> Result<[u8; N], Fail> {
        let bytes = self
            .data
            .get(self.pos..)
            .and_then(|rest| rest.get(..N))
            .and_then(|bytes| <[u8; N]>::try_from(bytes).ok())
            .ok_or(refuse(Rule::ReadPastEnd))?;
        self.pos += N;
        Ok(bytes)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:96 (`get_8`)
    fn u8(&mut self) -> Result<u8, Fail> {
        let [byte] = self.take::<1>()?;
        Ok(byte)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:925-937, 1164-1171
    fn u16(&mut self) -> Result<u16, Fail> {
        let bytes = self.take::<2>()?;
        Ok(if self.beq {
            u16::from_be_bytes(bytes)
        } else {
            u16::from_le_bytes(bytes)
        })
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:943-949, 1173-1180 (both
    // little-endian)
    fn u24(&mut self) -> Result<u32, Fail> {
        let [b0, b1, b2] = self.take::<3>()?;
        Ok(u32::from(b0) | u32::from(b1) << 8 | u32::from(b2) << 16)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:955-967, 1182-1189
    fn u32(&mut self) -> Result<u32, Fail> {
        let bytes = self.take::<4>()?;
        Ok(if self.beq {
            u32::from_be_bytes(bytes)
        } else {
            u32::from_le_bytes(bytes)
        })
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:973-985, 1191-1203 (both
    // little-endian)
    fn u48(&mut self) -> Result<u64, Fail> {
        let [b0, b1, b2, b3, b4, b5] = self.take::<6>()?;
        Ok(u64::from_le_bytes([b0, b1, b2, b3, b4, b5, 0, 0]))
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:991-1012, 1205-1217
    // ```c++
    //     bm::id64_t a = buf_[0]+
    //                    ((bm::id64_t)buf_[1] << 56) +
    //                    ((bm::id64_t)buf_[2] << 48) +
    //     ...
    //                    ((bm::id64_t)buf_[7] << 8);
    // ```
    fn u64(&mut self) -> Result<u64, Fail> {
        let bytes = self.take::<8>()?;
        if !self.beq {
            return Ok(u64::from_le_bytes(bytes));
        }
        let [b0, b1, b2, b3, b4, b5, b6, b7] = bytes;
        Ok(u64::from_be_bytes([b1, b2, b3, b4, b5, b6, b7, 0]) | u64::from(b0))
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:1103-1130 (`get_16(s, count)`:
    // a `memcpy` of `count` values), :1276-1295 (big-endian: a do-while loop that writes at least
    // one value, each from `buf_[0]` and `buf_[1]`, and moves 2 bytes for each)
    /// `get_16(out, count)`.
    fn u16s(&mut self, out: &mut [u16], count: u32) -> Result<(), Fail> {
        let count = count as usize;
        if self.beq {
            let n = count.max(1);
            let value = u16::from_be_bytes(self.take::<2>()?);
            out.get_mut(..n).ok_or(refuse(Rule::Bounds))?.fill(value);
            self.pos += 2 * (n - 1);
            return Ok(());
        }
        let bytes = self
            .data
            .get(self.pos..)
            .and_then(|rest| rest.get(..2 * count))
            .ok_or(refuse(Rule::ReadPastEnd))?;
        let out = out.get_mut(..count).ok_or(refuse(Rule::Bounds))?;
        for (value, pair) in out.iter_mut().zip(bytes.chunks_exact(2)) {
            *value = u16::from(pair[0]) | u16::from(pair[1]) << 8;
        }
        self.pos += 2 * count;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:1003-1027 (`get_32(w, count)`),
    // :1219-1236 (big-endian: a do-while loop that writes at least one word; the bit interval,
    // the one caller whose count may be 0, passes that count)
    /// `get_32(block + first, count)`, ORed into the block.
    fn words(&mut self, block: &mut Block, first: usize, count: usize) -> Result<(), Fail> {
        let bytes = self
            .data
            .get(self.pos..)
            .and_then(|rest| rest.get(..4 * count))
            .ok_or(refuse(Rule::ReadPastEnd))?;
        for (k, word) in bytes.chunks_exact(4).enumerate() {
            let word = [word[0], word[1], word[2], word[3]];
            let word = if self.beq {
                u32::from_be_bytes(word)
            } else {
                u32::from_le_bytes(word)
            };
            or_word32(block, first + k, word);
        }
        self.pos += 4 * count;
        Ok(())
    }
}

// ------------------------------------------------------------------------------------------------
// The bit reader.

/// `bm::bit_in` (encoding.h:298-409): a reader of bits, lowest first, that takes one 32-bit word
/// from the byte decoder each time it needs more bits; a new reader for each section, so the
/// bytes after a section start at its next word.
struct BitIn<'r, 'a> {
    src: &'r mut Bytes<'a>,
    acc: u32,
    used: u32,
}

/// `block_set_table<true>::_left[k]` (bmconst.h:371-376): the `k + 1` lowest bits.
fn left(k: u32) -> u32 {
    if k >= 31 {
        u32::MAX
    } else {
        (2 << k) - 1
    }
}

impl<'r, 'a> BitIn<'r, 'a> {
    fn new(src: &'r mut Bytes<'a>) -> Self {
        BitIn {
            src,
            acc: 0,
            used: 32,
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2638-2670 (`get_bits`; the
    // callers ask for 1 to 31 bits)
    fn get_bits(&mut self, count: u32) -> Result<u32, Fail> {
        debug_assert!((1..32).contains(&count));
        let mut acc = self.acc;
        let mut used = self.used;
        let free_bits = 32 - used;
        let value;
        if count <= free_bits || used == 32 {
            if count > free_bits {
                acc = self.src.u32()?;
                used = 0;
            }
            value = acc & (u32::MAX >> (32 - count));
            acc >>= count;
            used += count;
        } else {
            let low = acc;
            acc = self.src.u32()?;
            used = count - free_bits;
            value = low | (acc & (u32::MAX >> (32 - used))) << free_bits;
            acc >>= used;
        }
        self.acc = acc;
        self.used = used;
        Ok(value)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2675-2692 (`get_bit`)
    fn get_bit(&mut self) -> Result<u32, Fail> {
        if self.used != 32 {
            let value = self.acc & 1;
            self.acc >>= 1;
            self.used += 1;
            return Ok(value);
        }
        let word = self.src.u32()?;
        self.used = 1;
        self.acc = word >> 1;
        Ok(word & 1)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2506-2586 (`gamma`: the zero
    // bits, the border bit, then that many bits; `_left[zero_bits]` and `1 << zero_bits` with 32
    // or more zero bits are undefined, U3)
    fn gamma(&mut self) -> Result<u32, Fail> {
        let mut acc = self.acc;
        let mut used = self.used;
        if used == 32 {
            acc = self.src.u32()?;
            used = 0;
        }
        let mut zero_bits: u32 = 0;
        loop {
            if acc == 0 {
                zero_bits += 32 - used;
                if zero_bits > 31 {
                    return Err(refuse(Rule::GammaZeros));
                }
                used = 0;
                acc = self.src.u32()?;
                continue;
            }
            let first = acc.trailing_zeros();
            acc >>= first;
            zero_bits += first;
            used += first;
            break;
        }
        if zero_bits > 31 {
            return Err(refuse(Rule::GammaZeros));
        }
        // the border bit
        if used == 32 {
            acc = self.src.u32()?;
            used = 1;
        } else {
            used += 1;
        }
        acc >>= 1;
        let free_bits = 32 - used;
        let current;
        if zero_bits <= free_bits || used == 32 {
            if zero_bits > free_bits {
                acc = self.src.u32()?;
                used = 0;
            }
            current = (acc & left(zero_bits)) | 1 << zero_bits;
            acc >>= zero_bits;
            used += zero_bits;
        } else {
            let low = acc;
            acc = self.src.u32()?;
            used = zero_bits - free_bits;
            current = low | (acc & left(used)) << free_bits | 1 << zero_bits;
            acc >>= used;
        }
        self.acc = acc;
        self.used = used;
        Ok(current)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2441-2458 (`gamma8`)
    fn gamma8(&mut self) -> Result<u32, Fail> {
        Ok(match self.gamma()? {
            1 => self.gamma()?,
            2 => self.get_bits(8)?,
            3 => self.delta16()?,
            _ => 0,
        })
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2463-2483 (`delta16`)
    fn delta16(&mut self) -> Result<u32, Fail> {
        Ok(match self.gamma()? {
            1 => 511 - self.get_bits(8)?,
            2 => 512 + 255 - self.get_bits(8)?,
            3 => 512 + 256 + 255 - self.get_bits(8)?,
            _ => self.get_16_no()?,
        })
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2488-2501 (`delta16s`)
    fn delta16s(&mut self) -> Result<u32, Fail> {
        if self.get_bit()? == 1 {
            self.delta16()
        } else {
            self.get_bits(8)
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2591-2623 (`get_16_no`,
    // `get_24_no`, `get_32_no`: bytes of 8 bits, lowest first)
    fn get_16_no(&mut self) -> Result<u32, Fail> {
        let lo = self.get_bits(8)?;
        Ok(lo | self.get_bits(8)? << 8)
    }

    fn get_24_no(&mut self) -> Result<u32, Fail> {
        let v0 = self.get_bits(8)?;
        let v1 = self.get_bits(8)?;
        Ok(v0 + (v1 << 8) + (self.get_bits(8)? << 16))
    }

    fn get_32_no(&mut self) -> Result<u32, Fail> {
        let lo = self.get_16_no()?;
        Ok(lo | self.get_16_no()? << 16)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2177-2190 (the same lines in
    // each `bic_decode_*_cm*`)
    // ```c++
    //         unsigned val = hi - lo - sz + 1;
    //         if (val)
    //         {
    //             unsigned logv = bm::bit_scan_reverse32(val+1);
    //
    //             unsigned c = unsigned((1ull << (logv + 1)) - val - 1);
    //             int64_t half_c = c >> 1;
    //             int64_t half_r = val >> 1;
    //             int64_t lo1 = half_r - half_c - ((val + 1) & 1);
    //             int64_t hi1 = half_r + half_c + 1;
    //
    //             val = get_bits(logv);
    //             if (val <= lo1 || val >= hi1)
    //                 val += (get_bit() << logv);
    //         }
    // ```
    /// The centered minimal value of a range of `r + 1` values (always at most `r`; `r + 1 == 0`
    /// asks `bit_scan_reverse32(0)`, undefined).
    fn bic_value(&mut self, r: u32) -> Result<u32, Fail> {
        if r == 0 {
            return Ok(0);
        }
        let n = r.wrapping_add(1);
        if n == 0 {
            return Err(refuse(Rule::BicScan));
        }
        let logv = 31 - n.leading_zeros();
        let c = ((1u64 << (logv + 1)) - u64::from(r) - 1) as u32;
        let half_c = i64::from(c >> 1);
        let half_r = i64::from(r >> 1);
        let lo1 = half_r - half_c - i64::from(n & 1);
        let hi1 = half_r + half_c + 1;
        let mut x = self.get_bits(logv)?;
        if i64::from(x) <= lo1 || i64::from(x) >= hi1 {
            x += self.get_bit()? << logv;
        }
        Ok(x)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:336-341 (`bic_decode_u16`:
    // nothing for `sz` 0), :2209-2247 (`bic_decode_u16_cm`)
    // ```c++
    //         unsigned mid_idx = sz >> 1;
    //         val += lo + mid_idx;
    //         arr[mid_idx] = bm::gap_word_t(val);
    //
    //         if (sz <= 1)
    //             return;
    //
    //         bic_decode_u16_cm(arr, mid_idx, lo, bm::gap_word_t(val-1));
    //         arr += ++mid_idx;
    //         sz -= mid_idx;
    //         lo = bm::gap_word_t(val + 1);
    //     } while (sz);// for sz
    // ```
    /// `bic_decode_u16(arr, sz, lo, hi)`: the decoder writes `arr[0..sz]` (an array of NCBI's
    /// that is shorter is written past, U4).
    fn bic_u16(&mut self, arr: &mut [u16], sz: u32, lo: u16, hi: u16) -> Result<(), Fail> {
        if sz == 0 {
            return Ok(());
        }
        let arr = arr.get_mut(..sz as usize).ok_or(refuse(Rule::Bounds))?;
        self.bic_u16_cm(arr, lo, hi)
    }

    fn bic_u16_cm(&mut self, mut arr: &mut [u16], mut lo: u16, hi: u16) -> Result<(), Fail> {
        loop {
            let sz = arr.len() as u32;
            let r = u32::from(hi)
                .wrapping_sub(u32::from(lo))
                .wrapping_sub(sz)
                .wrapping_add(1);
            let val = self
                .bic_value(r)?
                .wrapping_add(u32::from(lo))
                .wrapping_add(sz >> 1);
            let (head, rest) = std::mem::take(&mut arr).split_at_mut((sz >> 1) as usize);
            let Some((mid, tail)) = rest.split_first_mut() else {
                return Ok(());
            };
            *mid = val as u16;
            if sz <= 1 {
                return Ok(());
            }
            self.bic_u16_cm(head, lo, val.wrapping_sub(1) as u16)?;
            arr = tail;
            lo = val.wrapping_add(1) as u16;
            if arr.is_empty() {
                return Ok(());
            }
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:343-348
    // (`bic_decode_u16_bitset`), :2251-2290 (`bic_decode_u16_cm_bitset`)
    // ```c++
    //         unsigned mid_idx = sz >> 1;
    //         val += lo + mid_idx;
    //
    //         // set bit in the target block
    //         {
    //             unsigned nword = (val >> bm::set_word_shift);
    //             block[nword] |= (1u << (val & bm::set_word_mask));
    //         }
    // ```
    /// `bic_decode_u16_bitset(block, sz, lo, hi)`: each value sets its bit (a value past the
    /// block's 65536 bits writes past it, U4); `dry`: `bic_decode_u16_dry` (encoding.h:352-357,
    /// 2306-2345), which NCBI runs for a full block: the same reads, no write.
    fn bic_u16_bitset(
        &mut self,
        block: &mut Block,
        sz: u32,
        lo: u16,
        hi: u16,
        dry: bool,
    ) -> Result<(), Fail> {
        if sz == 0 {
            return Ok(());
        }
        self.bic_u16_cm_bitset(block, sz, lo, hi, dry)
    }

    /// A range of `sz` values with `hi - lo + 1 == sz` holds `lo..=hi` and reads no bit, nor do its
    /// halves (each is such a range), so it is set at once: the decoder's time stays bounded by
    /// the bits it reads when `sz` is near 2^32 (a length below 2 that wraps).
    fn bic_u16_cm_bitset(
        &mut self,
        block: &mut Block,
        mut sz: u32,
        mut lo: u16,
        hi: u16,
        dry: bool,
    ) -> Result<(), Fail> {
        loop {
            let r = u32::from(hi)
                .wrapping_sub(u32::from(lo))
                .wrapping_sub(sz)
                .wrapping_add(1);
            if r == 0 && sz <= BLOCK_BITS {
                if !dry {
                    fill_range(block, u32::from(lo), u32::from(hi));
                }
                return Ok(());
            }
            let mid = sz >> 1;
            let val = self
                .bic_value(r)?
                .wrapping_add(u32::from(lo))
                .wrapping_add(mid);
            if !dry {
                if val >= BLOCK_BITS {
                    return Err(refuse(Rule::BicWrite));
                }
                set_bit(block, val);
            }
            if sz <= 1 {
                return Ok(());
            }
            self.bic_u16_cm_bitset(block, mid, lo, val.wrapping_sub(1) as u16, dry)?;
            sz -= mid + 1;
            lo = val.wrapping_add(1) as u16;
            if sz == 0 {
                return Ok(());
            }
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2167-2205
    // (`bic_decode_u32_cm`, called without the wrapper's check of `sz`: with `sz` 0 the loop's
    // body runs once and writes `arr[0]`)
    /// `bic_decode_u32_cm(arr, sz, lo, hi)`: writes `arr[0..max(sz, 1)]`.
    fn bic_u32(&mut self, arr: &mut [u32], sz: u32, lo: u32, hi: u32) -> Result<(), Fail> {
        if sz.max(1) as usize > arr.len() {
            return Err(refuse(Rule::Bounds));
        }
        self.bic_u32_cm(arr, sz, lo, hi)
    }

    fn bic_u32_cm(
        &mut self,
        mut arr: &mut [u32],
        mut sz: u32,
        mut lo: u32,
        hi: u32,
    ) -> Result<(), Fail> {
        loop {
            let r = hi.wrapping_sub(lo).wrapping_sub(sz).wrapping_add(1);
            let mid = sz >> 1;
            let val = self.bic_value(r)?.wrapping_add(lo).wrapping_add(mid);
            let (head, rest) = std::mem::take(&mut arr).split_at_mut(mid as usize);
            let Some((at, tail)) = rest.split_first_mut() else {
                return Err(refuse(Rule::Bounds));
            };
            *at = val;
            if sz <= 1 {
                return Ok(());
            }
            self.bic_u32_cm(head, mid, lo, val.wrapping_sub(1))?;
            arr = tail;
            sz -= mid + 1;
            lo = val.wrapping_add(1);
            if sz == 0 {
                return Ok(());
            }
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/encoding.h:2698-2798 (`decode_array`)
    /// `decode_array(arr, tb_wflags, &sz, default_sz)` into `arr` (`ex0_arr_` or `ex1_arr_`, of
    /// `EX_ENTRIES`): the header flags and the size. The forms write `arr[0..sz]` (delta-gamma
    /// also `arr[0]` with `sz` 0; a single value `arr[0]`); a window reads the array up to its
    /// width, which may pass `sz` into entries that an earlier token wrote.
    fn decode_array(
        &mut self,
        arr: &mut Arr,
        wflags: &mut Block,
        default_sz: u32,
    ) -> Result<(u32, u32), Fail> {
        let h3 = self.get_bits(8)?;
        if h3 & 3 == 3 && h3 & 0x80 != 0 {
            return Ok((h3, 0));
        }
        if h3 & 3 == 3 {
            let value = if h3 & 0x40 != 0 {
                0
            } else if h3 & 0x08 != 0 {
                self.gamma()? as u16
            } else {
                self.get_16_no()? as u16
            };
            arr.set(0, value)?;
            return Ok((h3, 1));
        }
        let sz = if default_sz != 0 {
            default_sz
        } else if h3 & 0x08 != 0 {
            self.gamma8()?.wrapping_add(1)
        } else {
            self.delta16()?
        };
        let min0 = if h3 & 0x40 != 0 {
            0
        } else {
            self.gamma()? as u16
        };
        let n = sz as usize;
        let a = arr.v.get_mut(..n.max(1)).ok_or(refuse(Rule::Bounds))?;
        if h3 & 3 == 0 {
            // delta-gamma
            let mut prev = if h3 & 0x80 != 0 {
                0
            } else {
                self.gamma()? as u16
            };
            a[0] = prev;
            for value in a.iter_mut().take(n).skip(1) {
                let delta = self.gamma()?;
                prev = u32::from(prev)
                    .wrapping_add(delta)
                    .wrapping_add(u32::from(min0)) as u16;
                *value = prev;
            }
            arr.wrote(n.max(1));
            return Ok((h3, sz));
        }
        if h3 & 2 != 0 {
            // gamma
            let zero_correct = u16::from(h3 & 0x80 != 0);
            for value in a.iter_mut().take(n) {
                *value = (self.gamma()? as u16)
                    .wrapping_sub(zero_correct)
                    .wrapping_add(min0);
            }
            arr.wrote(n);
            return Ok((h3, sz));
        }
        // BIC-DR
        if h3 & 0x80 != 0 {
            let mut min_v = self.get_16_no()? as u16;
            let mut max_v = self.get_16_no()? as u16;
            a[0] = min_v;
            *a.get_mut(n.wrapping_sub(1)).ok_or(refuse(Rule::Bounds))? = max_v;
            if sz != 2 {
                min_v = min_v.wrapping_add(1);
                max_v = max_v.wrapping_sub(1);
                self.bic_u16(&mut a[1..], sz.wrapping_sub(2), min_v, max_v)?;
            }
        } else {
            self.bic_u16(a, sz, 0, 65535)?;
        }
        arr.wrote(n);
        if sz == 2 && h3 & 0x80 != 0 {
            return Ok((h3, sz));
        }
        if self.get_bit()? != 0 {
            wflags.fill(0);
            let win_size = self.gamma()?;
            let wcnt = self.gamma()?.wrapping_add(15 - 1);
            let win_size = win_size.wrapping_add(9).wrapping_mul(2);
            // `*sz / win_size`: a division by zero (SIGFPE)
            if win_size == 0 {
                return Err(refuse(Rule::Window));
            }
            let max_wd = sz / win_size + 1;
            self.bic_u16_bitset(wflags, wcnt, 1, max_wd as u16, false)?;
            arr_restore_min_w(arr, sz, win_size, min0, wflags)?;
        } else if min0 != 0 {
            arr_restore_min_u16(&mut arr.v, sz, min0);
        }
        Ok((h3, sz))
    }
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:2648-2657 (`arr_restore_min`)
// ```c++
//     arr[0] = arr[0] + delta_acc;
//     for (unsigned i = 1; i < arr_len; ++i)
//     {
//         arr[i] += min0 + delta_acc;
//         BM_ASSERT(arr[i-1] < arr[i]);
//         delta_acc += min0;
//     } // for i
// ```
fn arr_restore_min_u16(arr: &mut [u16], len: u32, min0: u16) {
    let mut acc: u16 = 0;
    for value in arr.iter_mut().take(len as usize).skip(1) {
        *value = value.wrapping_add(min0).wrapping_add(acc);
        acc = acc.wrapping_add(min0);
    }
}

fn arr_restore_min_u32(arr: &mut [u32], min0: u32) {
    let mut acc: u32 = 0;
    for value in arr.iter_mut().skip(1) {
        *value = value.wrapping_add(min0).wrapping_add(acc);
        acc = acc.wrapping_add(min0);
    }
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:2517-2573 (`arr_restore_min_w`; its
// first loop reads and writes `arr[1..wlen]`, its second `arr[wlen-1..arr_len]`: past `arr_len`
// the entries are those that an earlier token wrote, and past the array NCBI reads past its buffer)
fn arr_restore_min_w(
    arr: &mut Arr,
    arr_len: u32,
    mut wlen: u32,
    min0: u16,
    wflags: &Block,
) -> Result<(), Fail> {
    let reach = arr_len.max(wlen) as usize;
    if reach > arr.v.len() {
        return Err(refuse(Rule::Bounds));
    }
    arr.check(reach)?;
    let arr = &mut arr.v[..reach];
    let mut delta_acc: u16 = 0;
    let mut min_w_prev: u32 = u32::MAX;
    for i in 1..wlen as usize {
        arr[i] = arr[i].wrapping_add(min0).wrapping_add(delta_acc);
        delta_acc = delta_acc.wrapping_add(min0);
        let delta = arr[i].wrapping_sub(arr[i - 1]);
        if u32::from(delta) < min_w_prev {
            min_w_prev = u32::from(delta);
        }
    }
    min_w_prev -= u32::from(min_w_prev != 0);
    let mut wave: u32 = 1;
    let mut i = wlen;
    while i < arr_len {
        if i + wlen > arr_len {
            wlen = arr_len % wlen;
        }
        let w_recalc = wflags[(wave >> 6) as usize & 1023] >> (wave & 63) & 1 != 0;
        let mut min_w: u32 = u32::MAX;
        for k in i as usize..(i + wlen) as usize {
            if w_recalc {
                let add = min_w_prev.wrapping_add(u32::from(delta_acc)) as u16;
                arr[k] = arr[k].wrapping_add(add);
                delta_acc = delta_acc.wrapping_add(min_w_prev as u16);
            } else {
                arr[k] = arr[k].wrapping_add(min0).wrapping_add(delta_acc);
                delta_acc = delta_acc.wrapping_add(min0);
            }
            let delta = arr[k].wrapping_sub(arr[k - 1]);
            if u32::from(delta) < min_w {
                min_w = u32::from(delta);
            }
        }
        min_w_prev = if min_w > u32::from(min0) {
            min_w - 1
        } else {
            u32::from(min0)
        };
        wave += 1;
        i += wlen;
    }
    Ok(())
}

// ------------------------------------------------------------------------------------------------
// Blocks.

fn set_bit(block: &mut Block, p: u32) {
    block[(p >> 6) as usize & 1023] |= 1 << (p & 63);
}

fn clear_bit(block: &mut Block, p: u32) {
    block[(p >> 6) as usize & 1023] &= !(1 << (p & 63));
}

fn or_word32(block: &mut Block, k: usize, word: u32) {
    block[(k >> 1) & 1023] |= u64::from(word) << ((k & 1) * 32);
}

/// Sets the bits `from..=to` (`to` below 65536).
fn fill_range(block: &mut Block, from: u32, to: u32) {
    let (first, last) = ((from >> 6) as usize, (to >> 6) as usize);
    let low = u64::MAX << (from & 63);
    let high = u64::MAX >> (63 - (to & 63));
    if first == last {
        block[first & 1023] |= low & high;
        return;
    }
    block[first & 1023] |= low;
    for word in block.iter_mut().take(last).skip(first + 1) {
        *word = u64::MAX;
    }
    block[last & 1023] |= high;
}

fn invert(block: &mut Block) {
    for word in block.iter_mut() {
        *word = !*word;
    }
}

/// The first bit at or after `from` that is set (`set`) or clear.
fn next_bit(block: &Block, from: u32, set: bool) -> Option<u32> {
    let flip = if set { 0 } else { u64::MAX };
    let mut w = (from >> 6) as usize;
    let mut word = (*block.get(w)? ^ flip) & (u64::MAX << (from & 63));
    loop {
        if word != 0 {
            return Some((w as u32) << 6 | word.trailing_zeros());
        }
        w += 1;
        word = *block.get(w)? ^ flip;
    }
}

/// The runs of set bits of `block`, or `None` above `limit` runs.
fn block_runs(block: &Block, limit: usize) -> Option<Vec<(u16, u16)>> {
    let mut runs = Vec::new();
    let mut from = 0;
    while let Some(first) = next_bit(block, from, true) {
        let end = next_bit(block, first, false).unwrap_or(BLOCK_BITS);
        if runs.len() == limit {
            return None;
        }
        runs.push((first as u16, (end - 1) as u16));
        if end == BLOCK_BITS {
            break;
        }
        from = end;
    }
    Some(runs)
}

// ------------------------------------------------------------------------------------------------
// GAP blocks.
//
// A GAP block is NCBI's array of `gap_word_t`: a header (`buf[0]`: the first run's value in bit 0,
// the level in bits 1-2, the last index `L` in bits 3-15), then the last bit of each run, the last
// of them 65535. LOSAT holds one as the entries that NCBI wrote and reads it entry by entry with
// NCBI's functions; an entry past those refuses the data (`Rule::Unwritten`). A GAP block "is one"
// when its ends increase up to 65535 (`gap_valid`); NCBI's searches of a block (`gap_bfind`,
// `gap_test_unr`, `gap_bit_count_range`, `gap_bit_count_unr`) run scalar, SSE2 or SSE4.2 code by
// the machine NCBI's BitMagic is built for, which agree on such a block. A query of one that is not
// runs all three (`gap_test3`, `gap_bfind3`, `gap_count3`) and is refused where they differ
// (`Rule::GapSearch`); so does a bit that the decoding sets in one (`gap_set_value`).

/// Entry `i` of a GAP block, of the entries that NCBI wrote.
fn gap_at(gap: &[u16], i: usize) -> Result<u16, Fail> {
    gap.get(i).copied().ok_or(refuse(Rule::Unwritten))
}

/// The GAP block's last index as its header holds it (`gap_length(buf) - 1`).
fn gap_last(gap: &[u16]) -> Result<usize, Fail> {
    Ok(usize::from(gap_at(gap, 0)? >> 3))
}

/// Whether the entries `gap` are a GAP block: a last index of at least 1 that the entries end at,
/// ends increasing up to 65535.
fn gap_valid(gap: &[u16]) -> bool {
    let Some((&header, ends)) = gap.split_first() else {
        return false;
    };
    usize::from(header >> 3) == ends.len()
        && ends.last() == Some(&65535)
        && ends.windows(2).all(|pair| pair[0] < pair[1])
}

/// The runs of 1s `(first, last)` of a GAP block that is one.
fn gap_runs(gap: &[u16]) -> Vec<(u16, u16)> {
    let mut runs = Vec::new();
    let mut value = gap[0] & 1 != 0;
    let mut from = 0u32;
    for &end in &gap[1..] {
        if value {
            runs.push((from as u16, end));
        }
        from = u32::from(end) + 1;
        value = !value;
    }
    runs
}

/// The GAP block of runs of 1s in order, apart and not touching: the form that NCBI's functions
/// keep a GAP block in (`gap_set_value`, `gap_init_range_block`), its header's level bits 0 and its
/// length the low 13 bits of the last index.
fn gap_of_runs(runs: &[(u16, u16)]) -> Vec<u16> {
    let mut gap = vec![u16::from(runs.first().is_some_and(|run| run.0 == 0))];
    for &(first, last) in runs {
        if first != 0 {
            gap.push(first - 1);
        }
        gap.push(last);
    }
    if gap.len() == 1 || gap[gap.len() - 1] != 65535 {
        gap.push(65535);
    }
    gap[0] |= ((gap.len() - 1) as u16).wrapping_shl(3);
    gap
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:4786-4811 (`gap_add_to_bitset`)
// ```c++
//     const T* pend = pcurr + len;
//     BM_ASSERT(*pend == 65535);
//     if (*pcurr & 1)  // Starts with 1
//     {
//         bm::or_bit_block(dest, 0, 1 + pcurr[1]);
//         pcurr += 3;
//     }
//     else
//         pcurr += 2;
//
//     unsigned bc, pos;
//     for (; pcurr <= pend; )
//     {
//         BM_ASSERT(*pcurr > pcurr[-1]);
//         pos = 1u + pcurr[-1];
//         bc = *pcurr - pcurr[-1];
//         pcurr += 2;
//         bm::or_bit_block(dest, pos, bc);
//     }
// ```
/// `gap_add_to_bitset(dest, gap, len)` (and `gap_convert_to_bitset`, which clears `dest` and
/// passes the header's length): ORs the runs of 1s whose end is at index `len` or before into the
/// block, with no other check on the GAP. A run that ends where it starts sets nothing
/// (`or_bit_block` with a count of 0 writes nothing, bmfunc.h:4517-4546). One that ends `d`
/// before it starts passes the count `2^32 - d`: `or_bit_block` then writes nothing when the
/// start's bit in its word (`bitpos`) is above `d` (the right margin `bitpos - d` is below 32 and
/// its mask meets none of the bits from `bitpos`); otherwise it shifts by 32 (`d == bitpos`) or
/// writes past the block: refused.
fn gap_add_to_bitset(gap: &[u16], len: usize, block: &mut Block) -> Result<(), Fail> {
    let mut k = 2;
    if gap_at(gap, 0)? & 1 != 0 {
        fill_range(block, 0, u32::from(gap_at(gap, 1)?));
        k = 3;
    }
    while k <= len {
        let (from, to) = (gap_at(gap, k - 1)?, gap_at(gap, k)?);
        if to < from {
            let bitpos = (u32::from(from) + 1) & 31;
            if u32::from(from - to) >= bitpos {
                return Err(refuse(Rule::GapOrder));
            }
            k += 2;
            continue;
        }
        if to > from {
            fill_range(block, u32::from(from) + 1, u32::from(to));
        }
        k += 2;
    }
    Ok(())
}

/// `gap_convert_to_bitset(block, gap)`: the bits of the GAP block by its header's length.
fn gap_bitmap(gap: &[u16]) -> Result<Box<Block>, Fail> {
    let mut block = Box::new([0; 1024]);
    gap_add_to_bitset(gap, gap_last(gap)?, &mut block)?;
    Ok(block)
}

/// The entries of a GAP block of each level (`gap_len_table`, bmconst.h:396-403).
const GAP_LEVELS: [usize; 4] = [128, 256, 512, 1280];

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:5409-5418 (`gap_calc_level` with
// the default `gap_len_table`: -1 above 1276)
/// `gap_calc_level(len, glevel_len)`: the level of a GAP block of `len` entries (`gap_length`,
/// one more than its last index).
fn gap_calc_level(len: usize) -> Option<u16> {
    GAP_LEVELS
        .iter()
        .position(|&n| len <= n - 4)
        .map(|level| level as u16)
}

/// Whether a GAP block of last index `len` is too long for a GAP block and converted to bits.
fn gap_converts(len: usize) -> bool {
    gap_calc_level(len + 1).is_none()
}

/// `gap_level(buf)` (bmfunc.h:1757-1761).
fn gap_level(head: u16) -> u16 {
    (head >> 1) & 3
}

/// `set_gap_level(buf, level)` (bmfunc.h:5387-5393).
fn with_level(head: u16, level: u16) -> u16 {
    (head & !6) | (level << 1)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:7414-7422 (`gap_operation_or`: the
// AND of both inverted, inverted), :3736-3790 (`gap_buff_op`)
// ```c++
//     T c1 = *cur1; T c2 = *cur2;
//     while (1)
//     {
//         bitval = (T) F::op(bitval1, bitval2);
//         res += (bitval != bitval_prev);
//         bitval_prev = bitval;
//         if (c1 < c2) // (*cur1 < *cur2)
//         {
//             *res = c1;
//             ++cur1; c1 = *cur1;
//             bitval1 ^= 1;
//         }
//         else // >=
//         {
//             *res = c2;
//             if (c2 < c1) // (*cur2 < *cur1)
//             {
//                 bitval2 ^= 1;
//             }
//             else  // equal
//             {
//                 if (c2 == (bm::gap_max_bits - 1))
//                     break;
//                 ++cur1; c1 = *cur1;
//                 bitval1 ^= 1; bitval2 ^= 1;
//             }
//             ++cur2; c2 = *cur2;
//         }
//     } // while
//     dlen = (unsigned)(res - dest);
//     *dest = (T)((*dest & 7) + (dlen << 3));
// ```
/// `gap_operation_or(a, b, tmp_buf, dlen)`: NCBI's merge of the two GAP blocks, which walks each to
/// the first place where both read 65535 (whatever their headers say), into the result buffer of
/// `OR_BUF` entries. The result's entries `[0..=dlen]`.
fn gap_or(a: &[u16], b: &[u16]) -> Result<Vec<u16>, Fail> {
    let mut bit1 = gap_at(a, 0)? & 1;
    let mut bit2 = gap_at(b, 0)? & 1;
    let mut prev = bit1 | bit2;
    let mut out = vec![prev, 0];
    let mut r = 1;
    let (mut i1, mut i2) = (1, 1);
    let mut c1 = gap_at(a, 1)?;
    let mut c2 = gap_at(b, 1)?;
    loop {
        let bit = bit1 | bit2;
        if bit != prev {
            r += 1;
            if r >= OR_BUF {
                return Err(refuse(Rule::Bounds));
            }
            out.push(0);
        }
        prev = bit;
        if c1 < c2 {
            out[r] = c1;
            i1 += 1;
            c1 = gap_at(a, i1)?;
            bit1 ^= 1;
        } else {
            out[r] = c2;
            if c2 < c1 {
                bit2 ^= 1;
            } else {
                if c2 == 65535 {
                    break;
                }
                i1 += 1;
                c1 = gap_at(a, i1)?;
                bit1 ^= 1;
                bit2 ^= 1;
            }
            i2 += 1;
            c2 = gap_at(b, i2)?;
        }
    }
    out[0] = (out[0] & 7).wrapping_add((r as u16).wrapping_shl(3));
    Ok(out)
}

/// `gap_is_all_zero` / `gap_is_all_one` (bmfunc.h): the first run, of 0s or of 1s, is the block.
fn gap_all(gap: &[u16], value: u16) -> bool {
    gap[0] & 1 == value && gap[1] == 65535
}

// ------------------------------------------------------------------------------------------------
// NCBI's searches of a GAP block that is not one.
//
// NCBI keeps such a block and reads it in the queries that reach it (`get_bit`, `count_range`).
// The search code is chosen when NCBI is built: scalar off x86-64, SSE2 on x86-64, SSE4.2 where
// NCBI's configure adds `-msse4.2` (configure.ac:2285-2330, ncbiconf_impl.h:326-338, bmdef.h:46-55,
// bmsimd.h:55-66). The three agree on a GAP block; on one that is not, LOSAT runs all three and
// answers where they agree (`Rule::GapSearch` where they differ: NCBI's answer depends on the
// machine it was built for). Each reads the block's entries with `gap_at` (an entry past the
// block's length is memory that NCBI did not write); a block of length 0 makes NCBI's binary
// searches run past the block or forever (`end - start` wraps), refused as `Rule::GapLength`.

/// The block's last index (at least 1).
fn gap_len(gap: &[u16]) -> Result<usize, Fail> {
    match gap_last(gap)? {
        0 => Err(refuse(Rule::GapLength)),
        len => Ok(len),
    }
}

/// One step of NCBI's binary searches: `if (buf[mid] < pos) start = mid + 1; else end = mid;`.
fn gap_step(gap: &[u16], pos: u32, start: &mut usize, end: &mut usize) -> Result<(), Fail> {
    let mid = (*start + *end) >> 1;
    if u32::from(gap_at(gap, mid)?) < pos {
        *start = mid + 1;
    } else {
        *end = mid;
    }
    Ok(())
}

/// The binary part that `gap_bfind`, `sse2_gap_bfind` and `sse42_gap_bfind` share: four steps
/// while 64 entries or more are left, then two while 16 or more (`step16` of them; one in
/// `sse42_gap_test`).
fn gap_narrow(
    gap: &[u16],
    pos: u32,
    mut start: usize,
    mut end: usize,
    step16: usize,
) -> Result<(usize, usize), Fail> {
    while end - start >= 64 {
        for _ in 0..4 {
            gap_step(gap, pos, &mut start, &mut end)?;
        }
    }
    while end - start >= 16 {
        for _ in 0..step16 {
            gap_step(gap, pos, &mut start, &mut end)?;
        }
    }
    Ok((start, end))
}

/// `is_set` of NCBI's searches: the value of the run that ends at entry `at`.
fn gap_value_at(gap: &[u16], at: usize) -> u16 {
    (gap[0] & 1) ^ ((at as u16).wrapping_sub(1) & 1)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:1899-1933 (`gap_test`, quoted in
// part)
// ```c++
//     unsigned start = 1;
//     unsigned end = 1 + ((*buf) >> 3);
//     if (end - start < 10)
//     {
//         unsigned sv = *buf & 1;
//         unsigned sv1= sv ^ 1;
//         if (buf[1] >= pos) return sv;
//         if (buf[2] >= pos) return sv1;
//         ...
//         if (buf[8] >= pos) return sv1;
//         BM_ASSERT(buf[9] >= pos);
//         return sv;
//     }
//     else
//     {
//         do
//         {
//             if (unsigned mid = (start + end) >> 1; buf[mid] < pos)
//                 start = mid + 1;
//             else
//                 end = mid;
//         } while (start != end);
//     }
//     return ((*buf) & 1) ^ ((--start) & 1);
// ```
/// `gap_test(buf, pos)`: `gap_test_unr` of NCBI's scalar build (bmfunc.h:1943-1955).
fn gap_test_scalar(gap: &[u16], pos: u32) -> Result<u16, Fail> {
    let len = gap_len(gap)?;
    if len < 10 {
        for k in 1..=8 {
            if u32::from(gap_at(gap, k)?) >= pos {
                return Ok(gap_value_at(gap, k));
            }
        }
        return Ok(gap[0] & 1);
    }
    let (mut start, mut end) = (1, 1 + len);
    while start != end {
        gap_step(gap, pos, &mut start, &mut end)?;
    }
    Ok(gap_value_at(gap, start))
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:1835-1886 (`gap_bfind` without
// `VECT_GAP_BFIND`: the binary part of `gap_narrow`, then)
// ```c++
//     for(; true; ++start)
//         if (buf[start] >= pos)
//             break;
//     *is_set = ((*buf) & 1) ^ ((start-1) & 1);
//     return start;
// ```
/// `gap_bfind(buf, pos, &is_set)` of NCBI's scalar build: the entry of the run that holds `pos`.
fn gap_bfind_scalar(gap: &[u16], pos: u32) -> Result<usize, Fail> {
    let len = gap_len(gap)?;
    let (mut start, _) = gap_narrow(gap, pos, 1, len, 2)?;
    while u32::from(gap_at(gap, start)?) < pos {
        start += 1;
    }
    Ok(start)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmsse2.h:1378-1430 (`sse2_gap_find`; the
// lanes from `size` on are set to 0xFFFF when `size` is below 8, none above), :1444-1499
// (`sse2_gap_bfind`: `end = 1 + len`, the binary part of `gap_narrow`,
// `size += (end != arr_end)`, `start += sse2_gap_find(buf + start, pos, size)`), :1501-1507
// (`sse2_gap_test`: its `is_set`)
/// `sse2_gap_bfind(buf, pos, &is_set)` of NCBI's SSE2 build.
fn gap_bfind_sse2(gap: &[u16], pos: u32) -> Result<usize, Fail> {
    let len = gap_len(gap)?;
    let arr_end = 1 + len;
    let (start, end) = gap_narrow(gap, pos, 1, arr_end, 2)?;
    let size = end - start + usize::from(end != arr_end);
    let at = |k: usize| gap_at(gap, start + k).map(u32::from);
    let found = if size < 4 {
        if at(0)? >= pos {
            0
        } else if at(1)? >= pos {
            1
        } else {
            2
        }
    } else {
        let mut found = None;
        for k in 0..8 {
            if (size < 8 && k >= size) || at(k)? >= pos {
                found = Some(k);
                break;
            }
        }
        match found {
            Some(k) => k,
            None if size == 8 => 8,
            None => {
                let base = size - 8;
                let mut k2 = size;
                for k in 0..8 {
                    if at(base + k)? >= pos {
                        k2 = base + k;
                        break;
                    }
                }
                k2
            }
        }
    };
    Ok(start + found)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmsse4.h:1339-1375 (`sse4_gap_find`: the
// count of lanes at or above `pos`, the lanes from `size` on set to 0xFFFF when `size` is below
// 8; `8 - count`, or with none, `size` less the count of the last 8 entries)
/// `sse4_gap_find(buf + start, pos, size)`, `size` 4 to 16.
fn sse4_gap_find(gap: &[u16], start: usize, pos: u32, size: usize) -> Result<usize, Fail> {
    let at = |k: usize| gap_at(gap, start + k).map(u32::from);
    let mut count = 0;
    for k in 0..8 {
        if (size < 8 && k >= size) || at(k)? >= pos {
            count += 1;
        }
    }
    if count != 0 {
        return Ok(8 - count);
    }
    let base = size - 8;
    let mut count = 0;
    for k in 0..8 {
        if at(base + k)? >= pos {
            count += 1;
        }
    }
    Ok(size - count)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmsse4.h:1389-1455 (`sse42_gap_bfind`),
// :1460-1526 (`sse42_gap_test`: one step per round while 16 or more are left)
// ```c++
//     unsigned start = 1;
//     unsigned end = ((*buf) >> 3);
//     unsigned size = end - start;
//     ... (the binary part)
//     ++size;
//     if (size < 4) // for very short vector use conventional scan
//     {
//         const unsigned short* BMRESTRICT pbuf = buf + start;
//         if (pbuf[0] >= pos) { }
//         else if (pbuf[1] >= pos) { start++; }
//         else
//         {
//             BM_ASSERT(pbuf[2] >= pos);
//             start+=2;
//         }
//     }
//     else
//     {
//         start += bm::sse4_gap_find(buf+start, (bm::gap_word_t)pos, size);
//     }
// ```
/// `sse42_gap_bfind` (`step16` 2) and `sse42_gap_test` (`step16` 1) of NCBI's SSE4.2 build.
fn gap_bfind_sse42(gap: &[u16], pos: u32, step16: usize) -> Result<usize, Fail> {
    let len = gap_len(gap)?;
    let (start, end) = gap_narrow(gap, pos, 1, len, step16)?;
    let size = end - start + 1;
    let at = |k: usize| gap_at(gap, start + k).map(u32::from);
    Ok(if size < 4 {
        if at(0)? >= pos {
            start
        } else if at(1)? >= pos {
            start + 1
        } else {
            start + 2
        }
    } else {
        start + sse4_gap_find(gap, start, pos, size)?
    })
}

/// The one answer of NCBI's three builds, or `Rule::GapSearch`.
fn agree<T: PartialEq>(answers: [T; 3]) -> Result<T, Fail> {
    let [a, b, c] = answers;
    if a == b && b == c {
        Ok(a)
    } else {
        Err(refuse(Rule::GapSearch))
    }
}

/// `gap_test_unr(buf, pos)` (bmfunc.h:1943-1955) of the three builds: `gap_test`,
/// `sse2_gap_test`, `sse42_gap_test`.
fn gap_test3(gap: &[u16], pos: u32) -> Result<bool, Fail> {
    agree([
        gap_test_scalar(gap, pos)?,
        gap_value_at(gap, gap_bfind_sse2(gap, pos)?),
        gap_value_at(gap, gap_bfind_sse42(gap, pos, 1)?),
    ])
    .map(|bit| bit != 0)
}

/// `gap_bfind(buf, pos, &is_set)` of the three builds (`VECT_GAP_BFIND`, bmfunc.h:1840-1841):
/// the entry of the run that holds `pos`.
fn gap_bfind3(gap: &[u16], pos: u32) -> Result<usize, Fail> {
    agree([
        gap_bfind_scalar(gap, pos)?,
        gap_bfind_sse2(gap, pos)?,
        gap_bfind_sse42(gap, pos, 2)?,
    ])
}

/// `*pcurr - prev` of NCBI's count functions: an `int` difference added to an `unsigned`.
fn gap_diff(to: u16, from: u32) -> u32 {
    u32::from(to).wrapping_sub(from)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:3097-3150 (`gap_bit_count_unr`; the
// scalar build's unrolled loop adds the same differences as its last loop), bmsse_util.h:1070-1091
// (`sse2_gap_sum_arr`, the SSE2 and SSE4.2 builds above 18 entries)
// ```c++
//     unsigned dsize = (*pcurr >> 3);
//     unsigned cnt = 0;
//     pcurr = buf + 1; // set up start position
//     T first_one = *buf & 1;
//     if (first_one)
//     {
//         cnt += *pcurr + 1;
//         ++pcurr;
//     }
//     ++pcurr;  // set GAP to 1
//     #elif defined(BMSSE42OPT) || defined(BMSSE2OPT)
//     if (dsize > 18)
//     {
//         const unsigned unr_factor = 16;
//         unsigned waves = (dsize - 2) / unr_factor;
//         pcurr = sse2_gap_sum_arr(pcurr, waves, &cnt);
//     }
//     ...
//     const T* pend = buf + dsize;
//     for ( ; pcurr <= pend ; pcurr+=2)
//         cnt += *pcurr - *(pcurr - 1);
// ```
/// `gap_bit_count_unr(buf)` of the scalar build (`sse` false) or of the SSE2 and SSE4.2 builds,
/// whose 16-bit lanes add the differences of 16 entries a wave modulo 2^16 in four sums.
fn gap_count_unr(gap: &[u16], sse: bool) -> Result<u32, Fail> {
    let len = gap_len(gap)?;
    let mut count = 0u32;
    let mut p = 1;
    if gap[0] & 1 != 0 {
        count = u32::from(gap_at(gap, 1)?) + 1;
        p = 2;
    }
    p += 1;
    if sse && len > 18 {
        let waves = (len - 2) / 16;
        let mut lanes = [0u16; 8];
        for _ in 0..waves {
            for (k, lane) in lanes.iter_mut().enumerate() {
                let sum = gap_at(gap, p + k - 1)?.wrapping_add(gap_at(gap, p + k + 7)?);
                *lane = lane.wrapping_add(sum);
            }
            p += 16;
        }
        for pair in lanes.chunks(2) {
            count = count.wrapping_add(u32::from(pair[1].wrapping_sub(pair[0])));
        }
    }
    while p <= len {
        count = count.wrapping_add(gap_diff(gap_at(gap, p)?, u32::from(gap_at(gap, p - 1)?)));
        p += 2;
    }
    Ok(count)
}

/// `block_bitcount` of a GAP block (bmblocks.h:1710-1720), the count of a whole block in
/// `count_range`, of the three builds.
fn gap_count3(gap: &[u16]) -> Result<u32, Fail> {
    let scalar = gap_count_unr(gap, false)?;
    let sse = gap_count_unr(gap, true)?;
    agree([scalar, sse, sse])
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:3163-3200 (`gap_bit_count_range`)
// ```c++
//     unsigned start_pos = bm::gap_bfind(buf, left, &is_set);
//     is_set = ~(is_set - 1u); // 0xFFF.. if true (mask for branchless code)
//     const T* pcurr = buf + start_pos;
//     if (right <= *pcurr) // we are in the target gap right now
//         bits_counter = unsigned(right - left + 1u) & is_set;
//     else
//     {
//         bits_counter = unsigned(*pcurr - left + 1u) & is_set;
//         ...
//             for (prev_gap = *pcurr++; right > *pcurr; prev_gap = *pcurr++)
//                 bits_counter += (is_set ^= ~0u) & (*pcurr - prev_gap);
//             bits_counter += unsigned(right - prev_gap) & (is_set ^ ~0u);
//     }
// ```
/// `gap_bit_count_range(buf, left, right)`, its `gap_bfind` of the three builds.
fn gap_count_range(gap: &[u16], left: u32, right: u32) -> Result<u32, Fail> {
    let mut p = gap_bfind3(gap, left)?;
    let mut mask = 0u32.wrapping_sub(u32::from(gap_value_at(gap, p)));
    let mut cur = u32::from(gap_at(gap, p)?);
    if right <= cur {
        return Ok(right.wrapping_sub(left).wrapping_add(1) & mask);
    }
    let mut count = cur.wrapping_sub(left).wrapping_add(1) & mask;
    let mut prev = cur;
    p += 1;
    cur = u32::from(gap_at(gap, p)?);
    while right > cur {
        mask = !mask;
        count = count.wrapping_add(mask & cur.wrapping_sub(prev));
        prev = cur;
        p += 1;
        cur = u32::from(gap_at(gap, p)?);
    }
    Ok(count.wrapping_add(right.wrapping_sub(prev) & !mask))
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:3488-3517 (`gap_bit_count_to`, the
// same in every build)
// ```c++
//     unsigned is_set = ~((unsigned(*buf) & 1u) - 1u); // 0xFFF.. if true (mask for branchless code)
//     const T* pcurr = buf + 1;
//     if (right <= *pcurr) // we are in the target block right now
//     {
//         bits_counter = (right + 1u) & is_set; // & is_set == if (is_set)
//     }
//     else
//     {
//         bits_counter = (*pcurr + 1u) & is_set;
//         prev_gap = *pcurr++;
//         for (is_set ^= ~0u; right > *pcurr; is_set ^= ~0u, prev_gap = *pcurr++)
//         {
//             bits_counter += (*pcurr - prev_gap) & is_set;
//             if (*pcurr == bm::gap_max_bits-1)
//                 goto cret;
//         }
//         bits_counter += (right - prev_gap) & is_set;
//     }
// ```
/// `gap_bit_count_to(buf, right)`.
fn gap_count_to(gap: &[u16], right: u32) -> Result<u32, Fail> {
    gap_len(gap)?;
    let mut mask = 0u32.wrapping_sub(u32::from(gap[0] & 1));
    let mut cur = u32::from(gap_at(gap, 1)?);
    if right <= cur {
        return Ok(right.wrapping_add(1) & mask);
    }
    let mut count = cur.wrapping_add(1) & mask;
    let mut prev = cur;
    let mut p = 2;
    mask = !mask;
    cur = u32::from(gap_at(gap, p)?);
    while right > cur {
        count = count.wrapping_add(cur.wrapping_sub(prev) & mask);
        mask = !mask;
        prev = cur;
        p += 1;
        cur = u32::from(gap_at(gap, p)?);
    }
    Ok(count.wrapping_add(right.wrapping_sub(prev) & mask))
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:4393-4438 (`gap_set_array`)
// ```c++
//     *buf = (T)((*buf & 6u) + (1u << 3)); // gap header setup
//     T* pcurr = buf + 1;
//     unsigned i = 0;
//     T curr = arr[i];
//     if (curr != 0) // need to add the first gap: (0 to arr[0]-1)
//     {
//         *pcurr = (T)(curr - 1);
//         ++pcurr;
//     }
//     else
//     {
//         ++(*buf); // GAP starts with 1
//     }
//     T prev = curr;
//     T acc = prev;
//     for (i = 1; i < len; ++i)
//     {
//         curr = arr[i];
//         if (curr == prev + 1)
//         {
//             ++acc;
//             prev = curr;
//         }
//         else
//         {
//             *pcurr++ = acc;
//             acc = curr;
//             *pcurr++ = (T)(curr-1);
//         }
//         prev = curr;
//     }
//     *pcurr = acc;
//     if (acc != bm::gap_max_bits - 1)
//     {
//         ++pcurr;
//         *pcurr = bm::gap_max_bits - 1;
//     }
//     unsigned gap_len = unsigned(pcurr - buf);
//     BM_ASSERT(gap_len == ((gap_len << 3) >> 3));
//     *buf = (T)((*buf & 7) + (gap_len << 3));
//     return gap_len+1;
// ```
/// `gap_set_array(gap, ids, len)` after `gap[0] = 0` (ids in any order; `arr[0]` is read also
/// for `len` 0, so `ids` holds at least one): the GAP's last index, whose header keeps its low 13
/// bits.
fn gap_set_array(gap: &mut Arr, ids: &[u16]) -> Result<usize, Fail> {
    fn put(gap: &mut [u16], at: usize, value: u16) -> Result<(), Fail> {
        *gap.get_mut(at).ok_or(refuse(Rule::Bounds))? = value;
        Ok(())
    }
    let (&first, rest) = ids.split_first().ok_or(refuse(Rule::Unwritten))?;
    let buf = &mut gap.v[..];
    buf[0] = 1 << 3;
    let mut at = 1;
    if first != 0 {
        put(buf, at, first - 1)?;
        at += 1;
    } else {
        buf[0] += 1;
    }
    let mut prev = first;
    let mut acc = prev;
    for &curr in rest {
        if u32::from(curr) == u32::from(prev) + 1 {
            acc = acc.wrapping_add(1);
        } else {
            put(buf, at, acc)?;
            acc = curr;
            put(buf, at + 1, curr.wrapping_sub(1))?;
            at += 2;
        }
        prev = curr;
    }
    put(buf, at, acc)?;
    if acc != 65535 {
        at += 1;
        put(buf, at, 65535)?;
    }
    buf[0] = (buf[0] & 7).wrapping_add((at as u16).wrapping_shl(3));
    gap.wrote(at + 1);
    Ok(at)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:4153-4205 (`gap_add_value`)
// ```c++
//     T end = (T)(*buf >> 3);
//     T curr = end;
//     T* pcurr = buf + end;
//     T* pend  = pcurr;
//     T* pprev = pcurr - 1;
//     if (!pos)
//     {
//         *buf ^= 1;
//         if ( buf[1] ) // We need to insert a 1 bit platform here.
//         {
//             ::memmove(&buf[2], &buf[1], (end - 1) * sizeof(gap_word_t));
//             buf[1] = 0;
//             ++end;
//         }
//         else // Only 1 bit in the GAP. We need to delete the first GAP.
//         {
//             pprev = buf + 1; pcurr = pprev + 1;
//             --end;
//             do { *pprev++ = *pcurr++; } while (pcurr < pend);
//         }
//     }
//     else if (((unsigned)(*pprev))+1 == pos && (curr > 1) ) // Left border bit
//     {
//        ++(*pprev);
//        if (*pprev == *pcurr)  // Curr. GAP to be merged with prev.GAP.
//             --end;
//     }
//     else if (*pcurr == pos) // Rightmost bit in the GAP. Border goes left.
//     {
//         --(*pcurr);
//         end += (pcurr == pend);
//     }
//     else  // Worst case we need to split current block.
//     {
//         pcurr[0] = (T)(pos-1);
//         pcurr[1] = (T)pos;
//         end = (T)(end+2);
//     }
//     *buf = (T)((*buf & 7) + (end << 3));
//     buf[end] = bm::gap_max_bits - 1;
// ```
/// `gap_add_value(gap, pos)`, step by step on `gap_temp_block_`: any order of ids, a header that
/// wraps (`end` is a `gap_word_t`). With `end` 0, the `memmove` of `end - 1` entries and the read
/// of `*pprev` before the array are undefined: refused.
fn gap_add_value(gap: &mut Arr, pos: u16) -> Result<(), Fail> {
    let mut end = gap.get(0)? >> 3;
    let cur = usize::from(end);
    if pos == 0 {
        gap.v[0] ^= 1;
        if gap.get(1)? != 0 {
            if cur == 0 {
                return Err(refuse(Rule::Bounds));
            }
            gap.copy_within(1..cur, 2)?;
            gap.set(1, 0)?;
            end += 1;
        } else {
            end = end.wrapping_sub(1);
            gap.copy_within(2..cur.max(3), 1)?;
        }
    } else {
        let prev = cur.checked_sub(1).ok_or(refuse(Rule::Bounds))?;
        let before = gap.get(prev)?;
        if u32::from(before) + 1 == u32::from(pos) && cur > 1 {
            let before = before.wrapping_add(1);
            gap.set(prev, before)?;
            if before == gap.get(cur)? {
                end -= 1;
            }
        } else if gap.get(cur)? == pos {
            gap.set(cur, pos.wrapping_sub(1))?;
            end += 1;
        } else {
            gap.set(cur, pos.wrapping_sub(1))?;
            gap.set(cur + 1, pos)?;
            end = end.wrapping_add(2);
        }
    }
    let header = (gap.v[0] & 7).wrapping_add(end.wrapping_shl(3));
    gap.set(0, header)?;
    gap.set(usize::from(end), 65535)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:4075-4140 (`gap_set_value`; the
// four-argument form of 4048-4058 runs `gap_set_value_cpos`, 3962-4038, the same steps)
// ```c++
//     unsigned is_set;
//     unsigned curr = bm::gap_bfind(buf, pos, &is_set);
//     T end = (T) (*buf >> 3);
//     if (is_set == val)
//         return end;
//     T* pcurr = buf + curr;
//     T* pprev = pcurr - 1;
//     T* pend = buf + end;
//     if (!pos)
//     {
//         *buf ^= 1;
//         if (buf[1]) // We need to insert a 1 bit GAP here
//         {
//             ::memmove(&buf[2], &buf[1], (end - 1) * sizeof(gap_word_t));
//             buf[1] = 0;
//             ++end;
//         }
//         else // Only 1 bit in the GAP. We need to delete the first GAP.
//         {
//             pprev = buf + 1; pcurr = pprev + 1;
//             goto copy_gaps;
//         }
//     }
//     else
//     if (curr > 1 && ((unsigned)(*pprev))+1 == pos) // Left border bit
//     {
//        ++(*pprev);
//        if (*pprev == *pcurr)  // Curr. GAP to be merged with prev.GAP.
//        {
//             --end;
//             if (pcurr != pend) // GAP merge: 2 GAPS to be deleted
//             {
//                 ++pcurr;
//                 copy_gaps:
//                 --end;
//                 do { *pprev++ = *pcurr++; } while (pcurr < pend);
//             }
//        }
//     }
//     else
//     if (*pcurr == pos) // Rightmost bit in the GAP. Border goes left.
//     {
//         --(*pcurr);
//         end += (pcurr == pend);
//     }
//     else  // Worst case: split current GAP
//     {
//         if (*pcurr != bm::gap_max_bits-1) // last gap does not need memmove
//             ::memmove(pcurr+2, pcurr, (end - curr + 1)*(sizeof(T)));
//         end += 2;
//         pcurr[0] = (T)(pos-1);
//         pcurr[1] = (T)pos;
//     }
//     *buf = (T)((*buf & 7) + (end << 3));
//     buf[end] = bm::gap_max_bits-1;
//     return end;
// ```
/// `gap_set_value(val, buf, pos)` after `gap_bfind` found the entry `curr`: any GAP block, a
/// header that wraps (`end` is a `gap_word_t`). A `memmove` whose count wraps (`end - 1` with
/// `end` 0, `end - curr + 1` with `curr` past `end + 1`) runs past the buffer: refused.
fn gap_set_at(gap: &mut Arr, val: u16, pos: u16, curr: usize) -> Result<u16, Fail> {
    let head = gap.get(0)?;
    let mut end = head >> 3;
    if gap_value_at(gap.written(), curr) == val {
        return Ok(end);
    }
    let pend = usize::from(end);
    if pos == 0 {
        gap.set(0, head ^ 1)?;
        if gap.get(1)? != 0 {
            if pend == 0 {
                return Err(refuse(Rule::Bounds));
            }
            gap.copy_within(1..pend, 2)?;
            gap.set(1, 0)?;
            end = end.wrapping_add(1);
        } else {
            end = end.wrapping_sub(1);
            gap.copy_within(2..pend.max(3), 1)?;
        }
    } else if curr > 1 && u32::from(gap.get(curr - 1)?) + 1 == u32::from(pos) {
        gap.set(curr - 1, pos)?;
        if pos == gap.get(curr)? {
            end = end.wrapping_sub(1);
            if curr != pend {
                end = end.wrapping_sub(1);
                gap.copy_within(curr + 1..pend.max(curr + 2), curr - 1)?;
            }
        }
    } else if gap.get(curr)? == pos {
        gap.set(curr, pos - 1)?;
        end = end.wrapping_add(u16::from(curr == pend));
    } else {
        if gap.get(curr)? != 65535 {
            if curr > pend + 1 {
                return Err(refuse(Rule::Bounds));
            }
            gap.copy_within(curr..pend + 1, curr + 2)?;
        }
        end = end.wrapping_add(2);
        gap.set(curr, pos - 1)?;
        gap.set(curr + 1, pos)?;
    }
    let head = (gap.get(0)? & 7).wrapping_add(end.wrapping_shl(3));
    gap.set(0, head)?;
    gap.set(usize::from(end), 65535)?;
    Ok(end)
}

/// `gap_set_value(val, buf, pos)` of NCBI's three builds, whose `gap_bfind` is `VECT_GAP_BFIND`
/// in the SIMD ones (bmfunc.h:1840-1841): the block as they leave it where they leave it the
/// same, else `Rule::GapSearch` (NCBI's block then depends on the machine it was built for). They
/// differ only on a GAP block that is not one.
fn gap_set_value(gap: &mut Arr, val: bool, pos: u16) -> Result<u16, Fail> {
    let found = {
        let entries = gap.written();
        let pos = u32::from(pos);
        [
            gap_bfind_scalar(entries, pos)?,
            gap_bfind_sse2(entries, pos)?,
            gap_bfind_sse42(entries, pos, 2)?,
        ]
    };
    if found[1..].iter().all(|&curr| curr == found[0]) {
        return gap_set_at(gap, u16::from(val), pos, found[0]);
    }
    // the bit already `val` in each build: none changes the block
    if found
        .iter()
        .all(|&curr| gap_value_at(gap.written(), curr) == u16::from(val))
    {
        return Ok(gap.get(0)? >> 3);
    }
    let mut builds = Vec::with_capacity(3);
    for curr in found {
        let mut build = gap.clone();
        let end = gap_set_at(&mut build, u16::from(val), pos, curr)?;
        builds.push((end, build));
    }
    let (end, first) = &builds[0];
    if builds[1..]
        .iter()
        .any(|(e, build)| e != end || build.written() != first.written())
    {
        return Err(refuse(Rule::GapSearch));
    }
    let end = *end;
    *gap = builds.swap_remove(0).1;
    Ok(end)
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:3000-3034 (`gap_restore_mins`)
fn gap_restore_mins(gap: &mut [u16], min0: u16, min1: u16) {
    let len = usize::from(gap[0] >> 3);
    gap[1] = gap[1].wrapping_add(min0);
    let mut acc = min0;
    for k in 2..len {
        if k % 2 == 0 {
            gap[k] = gap[k].wrapping_add(min1).wrapping_add(acc);
            acc = acc.wrapping_add(min1);
        } else {
            gap[k] = gap[k].wrapping_add(min0).wrapping_add(acc);
            acc = acc.wrapping_add(min0);
        }
    }
}

// ------------------------------------------------------------------------------------------------
// The vector being decoded.

/// A block of NCBI's vector while it decodes (`blocks_manager`, bmblocks.h); NCBI's null block is
/// no slot. The deserializer makes new blocks as GAP blocks (`set_new_blocks_strat(BM_GAP)`,
/// bmserial.h:5604).
enum Slot {
    /// A GAP block: the entries of its buffer that NCBI wrote (its header's level is the
    /// buffer's, `GAP_LEVELS`), and whether it is one (`gap_valid`: then its entries `[0..=L]`).
    Gap(Vec<u16>, bool),
    /// A GAP block that NCBI stored and cannot read back the same way each run (its header
    /// overwritten, `Decoder::gap_plain`). Refused when NCBI reads the block (a query that reaches
    /// it, `Kind::Poison`), freed unread when it makes the block full.
    Poison(Rule),
    /// A bitmap, also one of all ones (which is not the full block).
    Bits(Box<Block>),
    /// `FULL_BLOCK_FAKE_ADDR`, which NCBI neither reads nor writes.
    Full,
}

static FULL: Slot = Slot::Full;

impl Slot {
    /// A GAP block of the written entries `gap`: a block that is one keeps its entries
    /// `[0..=L]`, one that is not all that NCBI wrote (its searches may read past `L`).
    fn gap(mut gap: Vec<u16>) -> Slot {
        let len = gap.first().map_or(0, |&head| usize::from(head >> 3) + 1);
        if gap.get(..len).is_some_and(gap_valid) {
            gap.truncate(len);
            return Slot::Gap(gap, true);
        }
        Slot::Gap(gap, false)
    }
}

/// A sub-array of 256 blocks (`top_blocks_[i]`): null, `FULL_BLOCK_FAKE_ADDR` (256 full blocks),
/// or its blocks.
enum Lane {
    Absent,
    Full,
    Blocks(Box<[Option<Slot>; 256]>),
}

/// The blocks of NCBI's vector by sub-array, and how each token changes them: the operations
/// below are BitMagic's, named as NCBI's deserializer calls them.
struct Store {
    lanes: Vec<Lane>,
}

impl Default for Store {
    fn default() -> Self {
        Store {
            lanes: (0..256).map(|_| Lane::Absent).collect(),
        }
    }
}

fn or_into(block: &mut Block, bits: &Block) {
    for (word, bit) in block.iter_mut().zip(bits.iter()) {
        *word |= bit;
    }
}

/// `clone_gap_block(gap, gap_res)` (bmblocks.h:838-858): a new GAP block of the GAP's
/// `gap_length` entries (those that NCBI wrote; the rest of the copy is memory it did not write)
/// at the level of that length, or a bitmap above the last level.
fn clone_gap(gap: &[u16]) -> Result<Slot, Fail> {
    let len = gap_last(gap)?;
    let Some(level) = gap_calc_level(len + 1) else {
        return Ok(Slot::Bits(gap_bitmap(gap)?));
    };
    let mut block = gap[..gap.len().min(len + 1)].to_vec();
    block[0] = with_level(block[0], level);
    Ok(Slot::gap(block))
}

/// `assign_gap(res, res_len, blk)` (bmblocks.h:985-1030) into a GAP block of level `level`: a
/// bitmap above the last level (`convert_gap2bitset` with the result's length), a new block at
/// the result's level where it reaches the block's limit, else the block keeps its level.
fn assign_gap(mut res: Vec<u16>, level: u16) -> Result<Slot, Fail> {
    let len = res.len() - 1;
    let Some(new_level) = gap_calc_level(len + 1) else {
        let mut block = Box::new([0; 1024]);
        gap_add_to_bitset(&res, len, &mut block)?;
        return Ok(Slot::Bits(block));
    };
    let level = if len >= GAP_LEVELS[usize::from(level)] - 4 {
        new_level
    } else {
        level
    };
    res[0] = with_level(res[0], level);
    Ok(Slot::gap(res))
}

/// `extend_gap_block(nb, blk)` (bmblocks.h:1810-1836): a GAP block past its level's limit moves
/// to the next level (`allocate_gap_block` copies its `gap_length` entries and remakes the header),
/// or becomes a bitmap from the last level or at 1280 entries (`deoptimize_block`).
fn extend_gap(gap: &[u16]) -> Result<Slot, Fail> {
    let head = gap_at(gap, 0)?;
    let (level, len) = (gap_level(head), usize::from(head >> 3) + 1);
    if level == 3 || len >= GAP_LEVELS[3] {
        return Ok(Slot::Bits(gap_bitmap(gap)?));
    }
    let mut block = gap[..gap.len().min(len)].to_vec();
    block[0] = (((len - 1) as u16) << 3) | ((level + 1) << 1) | (head & 1);
    Ok(Slot::gap(block))
}

impl Store {
    fn get(&self, nb: u32) -> Option<&Slot> {
        match &self.lanes[(nb >> 8) as usize & 255] {
            Lane::Absent => None,
            Lane::Full => Some(&FULL),
            Lane::Blocks(blocks) => blocks[nb as usize & 255].as_ref(),
        }
    }

    /// `blk != 0`.
    fn exists(&self, nb: u32) -> bool {
        self.get(nb).is_some()
    }

    /// `IS_FULL_BLOCK(blk)`.
    fn is_full(&self, nb: u32) -> bool {
        matches!(self.get(nb), Some(Slot::Full))
    }

    /// The slot of block `nb`, its sub-array made (a full one as 256 full blocks).
    fn slot(&mut self, nb: u32) -> &mut Option<Slot> {
        let lane = &mut self.lanes[(nb >> 8) as usize & 255];
        if !matches!(lane, Lane::Blocks(_)) {
            let full = matches!(lane, Lane::Full);
            *lane = Lane::Blocks(Box::new(std::array::from_fn(|_| {
                full.then_some(Slot::Full)
            })));
        }
        match lane {
            Lane::Blocks(blocks) => &mut blocks[nb as usize & 255],
            _ => unreachable!(),
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bm.h:4646-4672 (`set_bit_no_check`),
    // bmblocks.h:1076-1118 (`check_allocate_block`: a null block becomes the empty GAP block of
    // `gap_set_all`), bm.h:2119-2132 (`gap_block_set_no_ret`: `gap_set_value`, then above the
    // level's limit `extend_gap_block`, bmblocks.h:1810-1836)
    /// `bvector<>::set_bit_no_check(nb * 65536 + p)`.
    fn set_bit(&mut self, nb: u32, p: u16) -> Result<(), Fail> {
        if self.is_full(nb) {
            return Ok(());
        }
        // `gap_set_all(allocate_gap_block(0), gap_max_bits, 0)`
        let slot = self
            .slot(nb)
            .get_or_insert_with(|| Slot::Gap(vec![1 << 3, 65535], true));
        match slot {
            Slot::Full => {}
            Slot::Poison(rule) => return Err(refuse(*rule)),
            Slot::Gap(gap, _) => {
                // `gap_block_set_no_ret`: `gap_set_value` in the block's buffer, then
                // `extend_gap_block` past the level's limit (`gap_limit`, bmfunc.h:1744-1748)
                let head = gap_at(gap, 0)?;
                let mut buf = Arr::block(gap, GAP_LEVELS[usize::from(gap_level(head))]);
                let end = usize::from(gap_set_value(&mut buf, true, p)?);
                let level = gap_level(buf.get(0)?);
                *slot = if usize::from(head >> 3) < end && end > GAP_LEVELS[usize::from(level)] - 4
                {
                    extend_gap(buf.written())?
                } else {
                    Slot::gap(buf.written().to_vec())
                };
            }
            Slot::Bits(block) => set_bit(block, u32::from(p)),
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmblocks.h:602-643
    // (`set_block_all_set`: the block is freed unread)
    fn set_block_all_set(&mut self, nb: u32) {
        if !matches!(self.lanes[(nb >> 8) as usize & 255], Lane::Full) {
            *self.slot(nb) = Some(Slot::Full);
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmblocks.h:644-702 (`set_all_set`: a
    // whole sub-array becomes `FULL_BLOCK_FAKE_ADDR` at once, or each of its blocks full)
    /// `set_all_set(nb, nb_to)`: blocks `nb..=nb_to` full.
    fn set_all_set(&mut self, nb: u32, nb_to: u32) {
        let (mut i_from, j_from) = (nb >> 8, nb & 255);
        let (mut i_to, j_to) = (nb_to >> 8, nb_to & 255);
        if i_from == i_to {
            self.set_lane_part(i_from, j_from, j_to);
            return;
        }
        if j_from > 0 {
            self.set_lane_part(i_from, j_from, 255);
            i_from += 1;
        }
        if j_to < 255 {
            self.set_lane_part(i_to, 0, j_to);
            i_to -= 1;
        }
        for i in i_from..=i_to {
            self.lanes[i as usize] = Lane::Full;
        }
    }

    /// Blocks `from..=to` of sub-array `i` full (`set_block_all_set_no_check` for each, unless
    /// the sub-array is full).
    fn set_lane_part(&mut self, i: u32, from: u32, to: u32) {
        if matches!(self.lanes[i as usize], Lane::Full) {
            return;
        }
        for j in from..=to {
            *self.slot(i << 8 | j) = Some(Slot::Full);
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bm.h:7835-7905 (`set_range_no_check`)
    // ```c++
    //     if ((nbit_left == 0) && (r == bm::bits_in_block - 1)) // full block
    //     {
    //         nb = nblock_left;
    //     }
    //     else
    //     {
    //         gap_init_range_block<gap_word_t>(tmp_gap_blk,
    //                                          (gap_word_t)nbit_left,
    //                                          (gap_word_t)r,
    //                                          (gap_word_t)1);
    //         ...
    //         combine_operation_with_block(nblock_left,
    //             BM_IS_GAP(block),
    //             block,
    //             (bm::word_t*) tmp_gap_blk,
    //             1, BM_OR);
    //
    //         if (nblock_left == nblock_right)  // in one block
    //             return;
    //         nb = nblock_left+1;
    //     }
    //     block_idx_type nb_to = nblock_right + (nbit_right ==(bm::bits_in_block-1));
    //     if (nb < nb_to)
    //     {
    //         blockman_.set_all_set(nb, nb_to-1);
    //     }
    //     if (nb_to > nblock_right)
    //         return;
    //     ...
    //     gap_init_range_block<gap_word_t>(tmp_gap_blk,
    //                                      (gap_word_t)0,
    //                                      (gap_word_t)nbit_right,
    //                                      (gap_word_t)1);
    //     combine_operation_with_block(nblock_right, ...);
    // ```
    /// `set_range_no_check(left, right)`, `left <= right`.
    fn set_range(&mut self, left: u32, right: u32) -> Result<(), Fail> {
        let (nb_left, nb_right) = (left >> 16, right >> 16);
        let nbit_left = (left & 0xFFFF) as u16;
        let nbit_right = (right & 0xFFFF) as u16;
        let r = if nb_left == nb_right {
            nbit_right
        } else {
            65535
        };
        let nb = if nbit_left == 0 && r == 65535 {
            nb_left
        } else {
            self.or_gap(nb_left, &gap_of_runs(&[(nbit_left, r)]), false)?;
            if nb_left == nb_right {
                return Ok(());
            }
            nb_left + 1
        };
        let nb_to = nb_right + u32::from(nbit_right == 65535);
        if nb < nb_to {
            self.set_all_set(nb, nb_to - 1);
        }
        if nb_to > nb_right {
            return Ok(());
        }
        self.or_gap(nb_right, &gap_of_runs(&[(0, nbit_right)]), false)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bm.h:7286-7360
    // (`combine_operation_block_or`), :7658-7833 (`combine_operation_with_block`, `BM_OR`)
    // ```c++
    //     if (IS_FULL_BLOCK(blk) || !arg_blk) // all bits are set
    //         return; // nothing to do
    //     ...
    //     if (BM_IS_GAP(blk)) // our block GAP-type
    //     {
    //         bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
    //         if (BM_IS_GAP(arg_blk))  // both blocks GAP-type
    //         {
    //             const bm::gap_word_t* res; unsigned res_len;
    //             res = bm::gap_operation_or(gap_blk, BMGAP_PTR(arg_blk),
    //                                        tmp_buf, res_len);
    //             blockman_.assign_gap_check(i, j, res, res_len, blk, tmp_buf);
    //             return;
    //         }
    //         // GAP or BIT
    //         //
    //         bm::word_t* new_blk = blockman_.get_allocator().alloc_bit_block();
    //         bm::bit_block_copy(new_blk, arg_blk);
    //         bm::gap_add_to_bitset(new_blk, gap_blk);
    //         ...
    //     }
    //     else // our block is BITSET-type (but NOT FULL_BLOCK we checked it)
    //     {
    //         if (BM_IS_GAP(arg_blk)) // argument block is GAP-type
    //         {
    //             const bm::gap_word_t* arg_gap = BMGAP_PTR(arg_blk);
    //             if (!blk)
    //             {
    //                 bool gap = true;
    //                 blk = blockman_.clone_gap_block(arg_gap, gap);
    //                 blockman_.set_block(i, j, blk, gap);
    //                 return;
    //             }
    //             // BIT & GAP
    //             bm::gap_add_to_bitset(blk, arg_gap);
    //             return;
    //         } // if arg_gap
    //     }
    //     ...
    //     bool all_one = bm::bit_block_or(blk, arg_blk);
    //     if (all_one)
    //     {
    //         BM_ASSERT(bm::is_bits_one((bm::wordop_t*) blk));
    //         blockman_.get_allocator().free_bit_block(blk);
    //         blockman_.set_block_ptr(i, j, FULL_BLOCK_FAKE_ADDR);
    //     }
    // ```
    /// ORs the GAP block `arg` (the entries NCBI wrote) into block `nb`: `block_or`
    /// (`combine_operation_block_or`, which makes a GAP result of all ones the full block), or
    /// `combine_operation_with_block` (whose only GAP arguments here are `set_range`'s ranges).
    fn or_gap(&mut self, nb: u32, arg: &[u16], block_or: bool) -> Result<(), Fail> {
        if self.is_full(nb) {
            return Ok(());
        }
        let slot = self.slot(nb);
        match slot {
            Some(Slot::Full) => {}
            None => *slot = Some(clone_gap(arg)?),
            Some(Slot::Poison(rule)) => return Err(refuse(*rule)),
            Some(Slot::Gap(gap, _)) => {
                // `assign_gap_check` (bmblocks.h:958-980), or `gap_is_all_zero` then `assign_gap`
                let level = gap_level(gap_at(gap, 0)?);
                let res = gap_or(gap, arg)?;
                *slot = if block_or && gap_all(&res, 1) {
                    Some(Slot::Full)
                } else if gap_all(&res, 0) {
                    None
                } else {
                    Some(assign_gap(res, level)?)
                };
            }
            Some(Slot::Bits(block)) => gap_add_to_bitset(arg, gap_last(arg)?, block)?,
        }
        Ok(())
    }

    /// ORs the bitmap `bits` into block `nb`: `block_or` (`combine_operation_block_or`, which
    /// makes a bitmap of all ones the full block), or `combine_operation_with_block`. A GAP block
    /// there becomes a bitmap (`gap_add_to_bitset` into a copy, `convert_gap2bitset`).
    fn or_bits(&mut self, nb: u32, bits: &Block, block_or: bool) -> Result<(), Fail> {
        if self.is_full(nb) {
            return Ok(());
        }
        let slot = self.slot(nb);
        match slot {
            Some(Slot::Full) => {}
            None => *slot = Some(Slot::Bits(Box::new(*bits))),
            Some(Slot::Poison(rule)) => return Err(refuse(*rule)),
            Some(Slot::Gap(gap, _)) => {
                let mut block = Box::new(*bits);
                gap_add_to_bitset(gap, gap_last(gap)?, &mut block)?;
                *slot = Some(Slot::Bits(block));
            }
            Some(Slot::Bits(block)) => {
                or_into(block, bits);
                if block_or && block.iter().all(|&word| word == u64::MAX) {
                    *slot = Some(Slot::Full);
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5326-5334, 5476-5492
    // (`decode_bit_block`, `decode_arrbit`), bmblocks.h:1574-1625 (`deoptimize_block`)
    /// The block that a bit form decodes into: a null block becomes an empty bitmap and a GAP
    /// block a bitmap; true for the full block, which stays (NCBI then decodes without writing).
    fn bit_target(&mut self, nb: u32) -> Result<bool, Fail> {
        if self.is_full(nb) {
            return Ok(true);
        }
        let slot = self.slot(nb);
        match slot {
            Some(Slot::Full) => return Ok(true),
            None => *slot = Some(Slot::Bits(Box::new([0; 1024]))),
            Some(Slot::Poison(rule)) => return Err(refuse(*rule)),
            Some(Slot::Gap(gap, _)) => *slot = Some(Slot::Bits(gap_bitmap(gap)?)),
            Some(Slot::Bits(_)) => {}
        }
        Ok(false)
    }

    /// ORs `bits` into the bitmap that `bit_target` made (`bit_block_or`, whose result is not
    /// checked for all ones); the full block first becomes a bitmap of ones (`deoptimize_block`,
    /// for the inverted forms).
    fn or_into_bits(&mut self, nb: u32, bits: &Block) {
        let slot = self.slot(nb);
        match slot {
            Some(Slot::Bits(block)) => or_into(block, bits),
            // the full block (`bit_target` made every other slot a bitmap)
            _ => *slot = Some(Slot::Bits(Box::new([u64::MAX; 1024]))),
        }
    }

    /// A new block `nb` (NCBI's `set_block` of a null block).
    fn put(&mut self, nb: u32, slot: Slot) {
        *self.slot(nb) = Some(slot);
    }

    /// The vector. A GAP block that is not one stays one, and a block whose entries the data does
    /// not decide a refusal, for the queries that reach them (`BitVector::get_bit`, `count_to`).
    fn finish(self, size: u32) -> Result<BitVector, Fail> {
        let mut segs = Segs::default();
        for (i, lane) in self.lanes.into_iter().enumerate() {
            let first = (i as u32) << 8;
            let blocks = match lane {
                Lane::Absent => continue,
                Lane::Full => {
                    segs.full(first, first + 255);
                    continue;
                }
                Lane::Blocks(blocks) => blocks,
            };
            for (j, slot) in blocks.into_iter().enumerate() {
                let nb = first | j as u32;
                match slot {
                    None => {}
                    Some(Slot::Full) => segs.full(nb, nb),
                    Some(Slot::Poison(rule)) => segs.poison(nb, rule),
                    Some(Slot::Gap(gap, false)) => segs.other_gap(nb, gap),
                    Some(Slot::Gap(gap, true)) => segs.runs(nb, gap_runs(&gap)),
                    Some(Slot::Bits(block)) => segs.bits(nb, block),
                }
            }
        }
        Ok(BitVector {
            size,
            segs: segs.segs,
            count: segs.count,
            other: segs.other,
            any_other: segs.any_other,
        })
    }
}

/// The decoded vector's segments as `finish` makes them: runs of full blocks, blocks of few runs,
/// blocks as words, GAP blocks that are not one.
struct Segs {
    segs: Vec<Seg>,
    count: u64,
    other: Result<u32, Fail>,
    any_other: bool,
}

impl Default for Segs {
    fn default() -> Self {
        Segs {
            segs: Vec::new(),
            count: 0,
            other: Ok(0),
            any_other: false,
        }
    }
}

impl Segs {
    fn push(&mut self, first: u32, last: u32, kind: Kind, bits: u64) {
        self.segs.push(Seg {
            first,
            last,
            before: self.count,
            other: self.other,
            kind,
        });
        self.count += bits;
    }

    /// A block that NCBI keeps as a GAP block that is not one, and its whole-block count.
    fn other_gap(&mut self, nb: u32, gap: Vec<u16>) {
        let count = gap_count3(&gap);
        self.push(nb, nb, Kind::Gap(gap.into_boxed_slice()), 0);
        self.add_other(count);
    }

    /// A block that a query of it refuses.
    fn poison(&mut self, nb: u32, rule: Rule) {
        self.push(nb, nb, Kind::Poison(rule), 0);
        self.add_other(Err(refuse(rule)));
    }

    fn add_other(&mut self, count: Result<u32, Fail>) {
        self.any_other = true;
        self.other = match (self.other, count) {
            (Ok(sum), Ok(count)) => Ok(sum.wrapping_add(count)),
            (Ok(_), Err(fail)) | (Err(fail), _) => Err(fail),
        };
    }

    fn full(&mut self, first: u32, last: u32) {
        let bits = u64::from(last - first + 1) * u64::from(BLOCK_BITS);
        if let Some(seg) = self.segs.last_mut() {
            if matches!(seg.kind, Kind::Full) && seg.last + 1 == first {
                seg.last = last;
                self.count += bits;
                return;
            }
        }
        self.push(first, last, Kind::Full, bits);
    }

    fn runs(&mut self, nb: u32, runs: Vec<(u16, u16)>) {
        if runs.is_empty() {
            return;
        }
        if runs == [(0, 65535)] {
            return self.full(nb, nb);
        }
        let mut before = Vec::with_capacity(runs.len() + 1);
        let mut sum = 0u32;
        before.push(sum);
        for &(a, b) in &runs {
            sum += u32::from(b) - u32::from(a) + 1;
            before.push(sum);
        }
        self.push(
            nb,
            nb,
            Kind::Runs(runs.into_boxed_slice(), before.into_boxed_slice()),
            u64::from(sum),
        );
    }

    fn bits(&mut self, nb: u32, words: Box<Block>) {
        if let Some(runs) = block_runs(&words, RUNS_MAX) {
            return self.runs(nb, runs);
        }
        let mut chunks = Box::new([0u32; 64]);
        let mut sum = 0u32;
        for (chunk, words) in chunks.iter_mut().zip(words.chunks(16)) {
            *chunk = sum;
            sum += words.iter().map(|word| word.count_ones()).sum::<u32>();
        }
        self.push(nb, nb, Kind::Bits(words, chunks), u64::from(sum));
    }
}

/// One of the deserializer's work arrays (bmserial.h:5110-5129), which NCBI allocates without
/// initializing them, once for each `bm::deserialize`, and the length of the prefix that this
/// decoding wrote (`hw`). An entry past it holds memory that NCBI did not write. Also a GAP block
/// while `gap_set_value` changes it (`Arr::block`).
#[derive(Clone)]
struct Arr {
    v: Vec<u16>,
    hw: usize,
}

impl Arr {
    fn new(len: usize) -> Arr {
        Arr {
            v: vec![0; len],
            hw: 0,
        }
    }

    /// A GAP block of `capacity` entries (`glevel_len[level]`, the allocation of
    /// `alloc_gap_block`) whose written entries are `gap`.
    fn block(gap: &[u16], capacity: usize) -> Arr {
        let mut v = gap[..gap.len().min(capacity)].to_vec();
        let hw = v.len();
        v.resize(capacity, 0);
        Arr { v, hw }
    }

    /// The entries `0..n` are written.
    fn wrote(&mut self, n: usize) {
        self.hw = self.hw.max(n.min(self.v.len()));
    }

    fn written(&self) -> &[u16] {
        &self.v[..self.hw]
    }

    /// Refused unless the entries `0..n` are written.
    fn check(&self, n: usize) -> Result<(), Fail> {
        if n > self.hw {
            return Err(refuse(Rule::Unwritten));
        }
        Ok(())
    }

    fn get(&self, i: usize) -> Result<u16, Fail> {
        gap_at(self.written(), i)
    }

    fn set(&mut self, i: usize, value: u16) -> Result<(), Fail> {
        *self.v.get_mut(i).ok_or(refuse(Rule::Bounds))? = value;
        if i == self.hw {
            self.hw += 1;
        }
        Ok(())
    }

    /// `memmove(&v[dest], &v[src.start], src.len())`: the copies of written entries are written,
    /// those of entries past the written prefix are not; a copy past the prefix leaves it.
    fn copy_within(&mut self, src: std::ops::Range<usize>, dest: usize) -> Result<(), Fail> {
        let len = src.len();
        if dest + len > self.v.len() || src.end > self.v.len() {
            return Err(refuse(Rule::Bounds));
        }
        if dest <= self.hw {
            let known = self.hw.saturating_sub(src.start).min(len);
            self.hw = if known == len {
                self.hw.max(dest + len)
            } else {
                dest + known
            };
        }
        self.v.copy_within(src, dest);
        Ok(())
    }

    /// `arr[0..max(n, 1)]` for a list of `n` ids (`gap_set_array` reads `arr[0]` also for 0).
    fn id_list(&self, n: usize) -> Result<&[u16], Fail> {
        let n = n.max(1);
        self.check(n)?;
        Ok(&self.v[..n])
    }
}

/// NCBI's work arrays of the deserializer (bmserial.h:5110-5129), and the block a token decodes
/// into (`temp_block_`).
struct Scratch {
    gap: Arr,
    ids: Arr,
    sb: Vec<u32>,
    ex0: Arr,
    ex1: Arr,
    wflags: Box<Block>,
    tmp: Box<Block>,
}

impl Default for Scratch {
    fn default() -> Self {
        Scratch {
            gap: Arr::new(GAP_ENTRIES),
            ids: Arr::new(GAP_ENTRIES),
            sb: vec![0; GAP_ENTRIES],
            ex0: Arr::new(EX_ENTRIES),
            ex1: Arr::new(EX_ENTRIES),
            wflags: Box::new([0; 1024]),
            tmp: Box::new([0; 1024]),
        }
    }
}

struct Decoder<'a> {
    bytes: Bytes<'a>,
    size: u32,
    store: Store,
    s: Scratch,
}

/// What the main loop does after a token: the next block (`++nb_i`), the same (`continue`), or
/// the end (`nb_i = bm::set_total_blocks`).
enum Go {
    Next,
    Stay,
    End,
}

// The serialization header's flags (bmserial.h:1157-1166).
const HM_RESIZE: u8 = 1 << 1;
const HM_ID_LIST: u8 = 1 << 2;
const HM_NO_BO: u8 = 1 << 3;
const HM_NO_GAPL: u8 = 1 << 4;
const HM_64_BIT: u8 = 1 << 5;

impl Decoder<'_> {
    /// `size_ = new_size` when a set bit or a range reaches past the size (`bvector<>::set_bit`,
    /// `set_range`: bm.h:2380-2384, 4232-4236).
    fn grow(&mut self, last: u32) {
        if last >= self.size {
            self.size = if last == u32::MAX { u32::MAX } else { last + 1 };
        }
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:4085-4125 (the free
    // `deserialize`: the byte-order byte picks the decoder), :5508-5991
    // (`deserializer::deserialize`), quoted in parts below.
    fn run(&mut self) -> Result<usize, Fail> {
        // ```c++
        //     bm::decoder dec(buf);
        //     unsigned char header_flag = dec.get_8();
        //     ByteOrder bo = bo_current;
        //     if (!(header_flag & BM_HM_NO_BO))
        //     {
        //         bo = (bm::ByteOrder) dec.get_8();
        //     }
        //     if (bo_current == bo)
        //     {
        //         deserializer<BV, bm::decoder> deserial;
        //     ...
        //     case LittleEndian:
        //         {
        //         deserializer<BV, bm::decoder_little_endian> deserial;
        // ```
        let flag = self.bytes.u8()?;
        if flag & HM_NO_BO == 0 {
            self.bytes.beq = self.bytes.u8()? != 1;
        }
        self.bytes.pos = 0;
        let flag = self.bytes.u8()?;
        if flag & HM_NO_BO == 0 {
            self.bytes.u8()?;
        }
        if flag & HM_64_BIT != 0 {
            return Err(Fail::Invalid);
        }
        // ```c++
        //     if (header_flag & BM_HM_ID_LIST)
        //     {
        //         // special case: the next comes plain list of integers
        //         if (header_flag & BM_HM_RESIZE)
        //         {
        //             ...
        //                 bv_size = dec.get_32();
        //             if (bv_size > bv.size())
        //                 bv.resize(bv_size);
        //         }
        //         for (unsigned cnt = dec.get_32(); cnt; --cnt)
        //         {
        //             bm::id_t idx = dec.get_32();
        //             bv.set(idx);
        //         } // for
        //         // -1 for compatibility with other deserialization branches
        //         return dec.size()-1;
        //     }
        // ```
        if flag & HM_ID_LIST != 0 {
            if flag & HM_RESIZE != 0 {
                let size = self.bytes.u32()?;
                self.size = self.size.max(size);
            }
            let count = self.bytes.u32()?;
            let mut ids = Vec::new();
            for _ in 0..count {
                let id = self.bytes.u32()?;
                self.grow(id);
                ids.push(id);
            }
            // `bv.set` comes before `set_new_blocks_strat(BM_GAP)` and makes bitmaps; the
            // decoding ends here, so only the bits count.
            ids.sort_unstable();
            for group in ids.chunk_by(|a, b| a >> 16 == b >> 16) {
                let tmp = &mut self.s.tmp;
                tmp.fill(0);
                for &id in group {
                    set_bit(tmp, id & (BLOCK_BITS - 1));
                }
                self.store.or_bits(group[0] >> 16, tmp, false)?;
            }
            return Ok(self.bytes.pos - 1);
        }
        if flag & HM_NO_GAPL == 0 {
            for _ in 0..4 {
                self.bytes.u16()?;
            }
        }
        if flag & HM_RESIZE != 0 {
            let size = self.bytes.u32()?;
            self.size = self.size.max(size);
        }
        let mut nb: u32 = 0;
        let mut bookmark: u32 = 0;
        loop {
            let token = self.bytes.u8()?;
            if token & 0x80 != 0 {
                nb = nb.wrapping_add(u32::from(token & 0x7F));
                if nb >= TOTAL_BLOCKS {
                    break;
                }
                continue;
            }
            let go = match token {
                0 | 9 => Go::End,
                1 => Go::Next,
                3 => {
                    nb = nb.wrapping_add(u32::from(self.bytes.u8()?));
                    Go::Stay
                }
                5 => {
                    nb = nb.wrapping_add(u32::from(self.bytes.u16()?));
                    Go::Stay
                }
                7 => {
                    nb = nb.wrapping_add(self.bytes.u32()?);
                    Go::Stay
                }
                10 => {
                    self.store.set_all_set(nb, TOTAL_BLOCKS - 1);
                    Go::End
                }
                2 => {
                    self.store.set_block_all_set(nb);
                    Go::Next
                }
                4 | 6 | 8 => {
                    let full_blocks = match token {
                        4 => u32::from(self.bytes.u8()?),
                        6 => u32::from(self.bytes.u16()?),
                        _ => self.bytes.u32()?,
                    };
                    self.full_blocks(nb, full_blocks)?;
                    nb = nb.wrapping_add(full_blocks.wrapping_sub(1));
                    Go::Next
                }
                26 => {
                    self.bytes.u64()?;
                    return Err(Fail::Invalid);
                }
                11 => {
                    // `decode_block_bit`: into a null block directly, else into `temp_block_`,
                    // ORed by `combine_operation_with_block`
                    let tmp = &mut self.s.tmp;
                    tmp.fill(0);
                    self.bytes.words(tmp, 0, 2048)?;
                    self.store.or_bits(nb, tmp, false)?;
                    Go::Next
                }
                19 => {
                    let bit = u32::from(self.bytes.u16()?).wrapping_add(nb << 16);
                    self.store.set_bit(bit >> 16, bit as u16)?;
                    Go::Next
                }
                22 => {
                    self.zero_runs(nb)?;
                    Go::Next
                }
                17 => {
                    self.interval(nb)?;
                    Go::Next
                }
                14 | 15 | 18 | 20 | 21 | 23 | 24 | 27 | 28 | 29 | 43 | 44 | 45 | 61 | 62 | 67 => {
                    self.gap_block(token, nb)?;
                    Go::Next
                }
                16 => {
                    self.arr_bit(nb)?;
                    Go::Next
                }
                30 | 31 | 32 | 33 | 34 | 57 | 63 | 64 | 65 | 66 => {
                    self.bit_block(token, nb)?;
                    Go::Next
                }
                56 | 68 => {
                    self.super_block(token)?;
                    nb = nb.wrapping_add(256 - (nb & 255));
                    Go::Stay
                }
                69 => {
                    nb = nb.wrapping_add(256 - (nb & 255));
                    Go::Stay
                }
                47..=49 => {
                    bookmark = nb;
                    match token {
                        47 => u32::from(self.bytes.u16()?),
                        48 => self.bytes.u24()?,
                        _ => self.bytes.u32()?,
                    };
                    Go::Stay
                }
                50..=55 => {
                    let sync = match token {
                        50 => u32::from(self.bytes.u8()?),
                        51 => u32::from(self.bytes.u16()?),
                        52 => self.bytes.u24()?,
                        53 => self.bytes.u32()?,
                        54 => self.bytes.u48()? as u32,
                        _ => self.bytes.u64()? as u32,
                    };
                    if nb != bookmark.wrapping_add(sync) {
                        return Err(Fail::Invalid);
                    }
                    Go::Stay
                }
                35..=42 | 58..=60 => return Err(refuse(Rule::Xor)),
                // 12, 13, 25, 46, 70-127
                _ => return Err(Fail::Invalid),
            };
            match go {
                Go::Next => nb = nb.wrapping_add(1),
                Go::Stay => {}
                Go::End => break,
            }
            if nb >= TOTAL_BLOCKS {
                break;
            }
        }
        Ok(self.bytes.pos)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5722-5729
    // ```c++
    //                 size_type from = nb_i * bm::gap_max_bits;
    //                 size_type to = from + full_blocks * bm::gap_max_bits;
    //                 bv.set_range(from, to-1);
    // ```
    // NCBI reference (598d8ae6): c++/include/util/bitset/bm.h:2368-2394 (`set_range`: the bounds in
    // order, the size grows to the right one)
    fn full_blocks(&mut self, nb: u32, full_blocks: u32) -> Result<(), Fail> {
        let from = nb.wrapping_mul(BLOCK_BITS);
        let to = from.wrapping_add(full_blocks.wrapping_mul(BLOCK_BITS));
        let (left, right) = (from, to.wrapping_sub(1));
        let (left, right) = if right < left {
            (right, left)
        } else {
            (left, right)
        };
        self.grow(right);
        self.store.set_range(left, right)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:4671-4698 (`read_0runs_block`)
    // ```c++
    //     unsigned char run_type = dec.get_8();
    //     for (unsigned j = 0; j < bm::set_block_size; run_type = !run_type)
    //     {
    //         unsigned run_length = dec.get_16();
    //         if (run_type)
    //         {
    //             unsigned run_end = j + run_length;
    //             BM_ASSERT(run_end <= bm::set_block_size);
    //             for (;j < run_end; ++j)
    //             {
    //                 unsigned w = dec.get_32();
    //                 blk[j] = w;
    //             }
    //         }
    //         else
    //         {
    //             j += run_length;
    //         }
    //     } // for j
    // ```
    fn zero_runs(&mut self, nb: u32) -> Result<(), Fail> {
        let tmp = &mut self.s.tmp;
        tmp.fill(0);
        let mut ones = self.bytes.u8()? != 0;
        let mut j: u32 = 0;
        while j < 2048 {
            let run = u32::from(self.bytes.u16()?);
            if ones {
                let end = j + run;
                if end > 2048 {
                    return Err(refuse(Rule::ZeroRuns));
                }
                self.bytes.words(tmp, j as usize, run as usize)?;
                j = end;
            } else {
                j += run;
            }
            ones = !ones;
        }
        // `combine_operation_with_block(nb_i, temp_block, 0, BM_OR)`
        self.store.or_bits(nb, tmp, false)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5443-5468
    // (`decode_block_bit_interval`)
    // ```c++
    //     unsigned head_idx = dec.get_16();
    //     unsigned tail_idx = dec.get_16();
    //     if (!blk)
    //     {
    //         ...
    //         for (unsigned k = 0; k < head_idx; ++k)
    //             blk[k] = 0;
    //         dec.get_32(blk + head_idx, tail_idx - head_idx + 1);
    //         for (unsigned j = tail_idx + 1; j < bm::set_block_size; ++j)
    //             blk[j] = 0;
    //     }
    //     else
    //     {
    //         bm::bit_block_set(temp_block_, 0);
    //         dec.get_32(temp_block_ + head_idx, tail_idx - head_idx + 1);
    //         bv.combine_operation_with_block(nb, temp_block_, 0, BM_OR);
    //     }
    // ```
    /// Token 17. `get_32` moves `tail - head + 1` words (wrapping), the big-endian decoder at
    /// least one; with `tail == head - 1` that one word goes to word `head`, which a null block's
    /// second loop clears again. Refused where the words or a null block's first loop pass the
    /// block's 2048 words.
    fn interval(&mut self, nb: u32) -> Result<(), Fail> {
        let head = u32::from(self.bytes.u16()?);
        let tail = u32::from(self.bytes.u16()?);
        let count = tail.wrapping_sub(head).wrapping_add(1);
        let words = if self.bytes.beq { count.max(1) } else { count };
        let exists = self.store.exists(nb);
        if (!exists && head > 2048) || (words > 0 && u64::from(head) + u64::from(words) > 2048) {
            return Err(refuse(Rule::Interval));
        }
        let tmp = &mut self.s.tmp;
        tmp.fill(0);
        self.bytes.words(tmp, head as usize, words as usize)?;
        if !exists && count == 0 {
            tmp.fill(0);
        }
        self.store.or_bits(nb, tmp, false)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5154-5317
    // (`deserialize_gap`), :4744-5034 (`read_gap_block`), :4200-4278 (`read_id_list`); quoted at
    // each form below. Each form makes a GAP block in `gap_temp_block_` (14 and 15 into a null
    // block: the new block itself), which `combine_operation_block_or` ORs into block `nb`; the id
    // lists of `deserialize_gap` above the last level are converted to bits first.
    // ```c++
    //     else
    //     {
    //         bm::word_t* tmp_blk = (bm::word_t*)gap_temp_block;
    //         BMSET_PTRGAP(tmp_blk);
    //         bv.combine_operation_block_or(i0, j0, blk, tmp_blk);
    //     }
    // ```
    fn gap_block(&mut self, token: u8, nb: u32) -> Result<(), Fail> {
        let s = &mut self.s;
        let gap = &mut s.gap;
        match token {
            14 | 15 => return self.gap_plain(nb),
            61 => return self.gap_v3(nb),
            18 | 21 | 23 | 28 | 29 | 44 | 45 => {
                // ```c++
                //         	unsigned arr_len = this->read_id_list(dec, btype, this->id_array_);
                //             gap_temp_block[0] = 0; // reset unused bits in gap header
                //             unsigned gap_len =
                //               bm::gap_set_array(gap_temp_block, this->id_array_, arr_len);
                //             if (inv_flag)
                //                 bm::gap_invert(gap_temp_block);
                //
                //             BM_ASSERT(gap_len == bm::gap_length(gap_temp_block));
                //             int level = bm::gap_calc_level(gap_len, bman.glen());
                //             if (level == -1)  // Too big to be GAP: convert to BIT block
                //             {
                //                 bm::gap_convert_to_bitset(temp_block_, gap_temp_block);
                //                 bv.combine_operation_block_or(i0, j0, blk, temp_block_);
                //                 return;
                //             }
                // ```
                // (23 and 29 in `read_gap_block`: the same, inverted, and no conversion.) The level
                // is of the GAP's length, the conversion of the header's 13 bits of it.
                let n = read_id_list(&mut self.bytes, token, &mut s.ids)?;
                let full = gap_set_array(gap, s.ids.id_list(n)?)?;
                if matches!(token, 23 | 29 | 45) {
                    gap.v[0] ^= 1;
                }
                if !matches!(token, 23 | 29) && gap_converts(full) {
                    let tmp = &mut s.tmp;
                    tmp.fill(0);
                    gap_add_to_bitset(gap.written(), usize::from(gap.v[0] >> 3), tmp)?;
                    return self.store.or_bits(nb, tmp, true);
                }
            }
            24 => {
                // ```c++
                //             bm::gap_set_all(dst_block, bm::gap_max_bits, 0);
                //             gap_word_t len = decoder.get_16();
                //             for (gap_word_t k = 0; k < len; ++k)
                //             {
                //                 gap_word_t bit_idx = decoder.get_16();
                // 				bm::gap_add_value(dst_block, bit_idx);
                //             } // for
                //     ...
                //         bm::gap_invert(dst_block);
                // ```
                // `gap_set_all` keeps the header's level bits, which no step here reads.
                gap.set(0, (gap.v[0] & 6) + (1 << 3))?;
                gap.set(1, 65535)?;
                let n = self.bytes.u16()?;
                for _ in 0..n {
                    let id = self.bytes.u16()?;
                    gap_add_value(gap, id)?;
                }
                gap.v[0] ^= 1;
            }
            20 => {
                // ```c++
                //         unsigned len = (gap_head >> 3);
                //         --len;
                //         *dst_block = gap_head;
                //         gap_word_t* gap_data_ptr = dst_block + 1;
                //         bit_in_type bin(decoder);
                //         gap_word_t v = (gap_word_t)bin.gamma();
                //         gap_word_t gap_sum = *gap_data_ptr = (gap_word_t)(v - 1);
                //         for (unsigned i = 1; i < len; ++i)
                //         {
                //             v = (gap_word_t)bin.gamma();
                //             gap_sum = (gap_word_t)(gap_sum + v);
                //             *(++gap_data_ptr) = gap_sum;
                //         }
                //         dst_block[len+1] = bm::gap_max_bits - 1;
                // ```
                // A length of 0 runs the loop 2^32 - 1 times, past `gap_temp_block_`.
                let head = self.bytes.u16()?;
                let len = usize::from(head >> 3);
                if len == 0 {
                    return Err(refuse(Rule::GapLength));
                }
                gap.v[0] = head;
                let mut bin = BitIn::new(&mut self.bytes);
                let mut sum = (bin.gamma()? as u16).wrapping_sub(1);
                gap.v[1] = sum;
                for k in 2..len {
                    sum = sum.wrapping_add(bin.gamma()? as u16);
                    gap.v[k] = sum;
                }
                gap.v[len] = 65535;
                gap.wrote(len + 1);
            }
            27 => {
                // ```c++
                //         unsigned len = (gap_head >> 3);
                //         *dst_block = gap_head;
                //         bm::gap_word_t min_v = decoder.get_16();
                //         dst_block[1] = min_v;
                //         bit_in_type bin(decoder);
                //         bin.bic_decode_u16(&dst_block[2], len-2, min_v, 65535);
                //         dst_block[len] = bm::gap_max_bits - 1;
                // ```
                let head = self.bytes.u16()?;
                let len = usize::from(head >> 3);
                gap.v[0] = head;
                let min = self.bytes.u16()?;
                gap.v[1] = min;
                BitIn::new(&mut self.bytes).bic_u16(
                    &mut gap.v[2..],
                    (len as u32).wrapping_sub(2),
                    min,
                    65535,
                )?;
                gap.v[len] = 65535;
                gap.wrote(len + 1);
            }
            43 => {
                // ```c++
                //             unsigned len = (gap_head >> 3);
                //             bm::gap_word_t min8 = gap_head & bm::h2f_min_v_8bit;
                //             bm::gap_word_t tail8 = gap_head & bm::h2f_max_v_8bit;
                //             gap_head &= bm::gap_word_t(~(3 << 1)); // clear the flags
                //             if (min8)
                //                 min_v = decoder.get_8();
                //             else
                //                 min_v = decoder.get_16();
                //             if (tail8)
                //                 max_v = decoder.get_8();
                //             else
                //                 max_v = decoder.get_16();
                //             max_v = bm::gap_word_t(65535 - max_v); // tail correction
                //             dst_block[0] = gap_head;
                //             dst_block[1] = min_v;
                //             bit_in_type bin(decoder);
                //             bin.bic_decode_u16(&dst_block[2], len-3, min_v, max_v);
                //             dst_block[len-1] = max_v;
                //             dst_block[len] = bm::gap_max_bits - 1;
                // ```
                let head = self.bytes.u16()?;
                let len = usize::from(head >> 3);
                let min = if head & 2 != 0 {
                    u16::from(self.bytes.u8()?)
                } else {
                    self.bytes.u16()?
                };
                let tail = if head & 4 != 0 {
                    u16::from(self.bytes.u8()?)
                } else {
                    self.bytes.u16()?
                };
                let max = 65535 - tail;
                gap.v[0] = head & !6;
                gap.v[1] = min;
                BitIn::new(&mut self.bytes).bic_u16(
                    &mut gap.v[2..],
                    (len as u32).wrapping_sub(3),
                    min,
                    max,
                )?;
                gap.v[len - 1] = max;
                gap.v[len] = 65535;
                gap.wrote(len + 1);
            }
            62 => {
                // ```c++
                //             gap_head = (bm::gap_word_t) bin.delta16s();
                //             unsigned len = (gap_head >> 3);
                //             ...
                //             min_v = (bm::gap_word_t) ((min8) ? bin.gamma8() : bin.get_16_no());
                //             max_v = (bm::gap_word_t) ((tail8) ? bin.gamma8() :  bin.get_16_no());
                //             max_v = bm::gap_word_t(65535 - max_v); // tail correction
                //             dst_block[0] = gap_head;
                //             dst_block[1] = min_v;
                //             bin.bic_decode_u16(&dst_block[2], len-3, min_v, max_v);
                //             dst_block[len-1] = max_v;
                //             dst_block[len] = bm::gap_max_bits - 1;
                // ```
                let mut bin = BitIn::new(&mut self.bytes);
                let head = bin.delta16s()? as u16;
                let len = usize::from(head >> 3);
                let min = if head & 2 != 0 {
                    bin.gamma8()?
                } else {
                    bin.get_16_no()?
                } as u16;
                let tail = if head & 4 != 0 {
                    bin.gamma8()?
                } else {
                    bin.get_16_no()?
                } as u16;
                let max = 65535 - tail;
                gap.v[0] = head & !6;
                gap.v[1] = min;
                bin.bic_u16(&mut gap.v[2..], (len as u32).wrapping_sub(3), min, max)?;
                gap.v[len - 1] = max;
                gap.v[len] = 65535;
                gap.wrote(len + 1);
            }
            _ => {
                // 67
                // ```c++
                //         unsigned len = bin.gamma() + 1;
                //         unsigned start_flag = bin.get_bit();
                //         unsigned use_gamma = bin.get_bit();
                //         dst_block[0] = gap_head = bm::gap_word_t((len << 3) | start_flag);
                //         if (use_gamma)
                //         {
                //             dst_block[1] = (bm::gap_word_t) bin.gamma8();
                //             bm::gap_word_t prev = dst_block[1];
                //             for( unsigned i = 2; i < len; ++i)
                //             {
                //                 unsigned delta = bin.gamma8();
                //                 dst_block[i] = (bm::gap_word_t)(prev + delta);
                //                 prev = dst_block[i];
                //             } // for i
                //         }
                //         else
                //         {
                //             for(unsigned i = 1; i < len; ++i)
                //             {
                //                 dst_block[i] = (bm::gap_word_t) bin.get_16_no();
                //             } // for i
                //         }
                //         dst_block[len] = bm::gap_max_bits - 1;
                // ```
                // A length of 0 (a gamma value of 2^32 - 1) writes 65535 into the header; one
                // above 65535 writes past `gap_temp_block_`. The header keeps the low 13 bits of
                // the length.
                let mut bin = BitIn::new(&mut self.bytes);
                let len = bin.gamma()?.wrapping_add(1) as usize;
                let start = bin.get_bit()?;
                let use_gamma = bin.get_bit()?;
                if len >= GAP_ENTRIES {
                    return Err(refuse(Rule::Bounds));
                }
                gap.set(0, ((len as u16).wrapping_shl(3)) | start as u16)?;
                if use_gamma != 0 {
                    let mut prev = bin.gamma8()? as u16;
                    gap.set(1, prev)?;
                    for k in 2..len {
                        prev = u32::from(prev).wrapping_add(bin.gamma8()?) as u16;
                        gap.v[k] = prev;
                    }
                } else {
                    for k in 1..len {
                        gap.v[k] = bin.get_16_no()? as u16;
                    }
                }
                gap.v[len] = 65535;
                gap.wrote(len + 1);
            }
        }
        self.store.or_gap(nb, self.s.gap.written(), true)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5174-5229 (`deserialize_gap`,
    // `set_block_gap`, `set_block_gapbit`)
    // ```c++
    //         gap_head = (gap_word_t)
    //             (sizeof(gap_word_t) == 2 ? dec.get_16() : dec.get_32());
    //         unsigned len = gap_length(&gap_head);
    //         int level = gap_calc_level(len, bman.glen());
    //         --len;
    //         if (level == -1)  // Too big to be GAP: convert to BIT block
    //         {
    //             *gap_temp_block = gap_head;
    //             dec.get_16(gap_temp_block+1, len - 1);
    //             gap_temp_block[len] = gap_max_bits - 1;
    //             if (blk == 0)  // block does not exist yet
    //             {
    //                 ...
    //                 gap_convert_to_bitset(blk, gap_temp_block);
    //             }
    //             else  // We have some data already here. Apply OR operation.
    //             {
    //                 gap_convert_to_bitset(temp_block_,
    //                                       gap_temp_block);
    //                 bv.combine_operation_block_or(i0, j0, blk, temp_block_);
    //             }
    //             return;
    //         } // level == -1
    //         ...
    //         if (blk == 0)
    //         {
    //             ...
    //             *gap_blk_ptr = gap_head;
    //             ...
    //             if (len > 1)
    //                 dec.get_16(gap_blk + 1, len - 1);
    //             gap_blk[len] = bm::gap_max_bits - 1;
    //         }
    //         else // target block exists
    //         {
    //             // read GAP block into a temp memory and perform OR
    //             *gap_temp_block = gap_head;
    //             dec.get_16(gap_temp_block + 1, len - 1);
    //             gap_temp_block[len] = bm::gap_max_bits - 1;
    //             break;
    //         }
    // ```
    /// Tokens 14 and 15. Into a null block the GAP becomes the block (a length of 0 puts 65535 in
    /// the header of a 128-entry block: NCBI cannot read it back). Into a block that exists (the
    /// full block too) it is read into `gap_temp_block_`: a length of 0 reads 2^32 - 1 values,
    /// refused; a length of 1 reads 2 bytes with the big-endian decoder and none with the
    /// little-endian one.
    fn gap_plain(&mut self, nb: u32) -> Result<(), Fail> {
        let head = self.bytes.u16()?;
        let len = usize::from(head >> 3);
        let exists = self.store.exists(nb);
        let gap = &mut self.s.gap;
        if gap_converts(len) {
            gap.v[0] = head;
            self.bytes.u16s(&mut gap.v[1..], (len - 1) as u32)?;
            gap.v[len] = 65535;
            gap.wrote(len + 1);
            let tmp = &mut self.s.tmp;
            tmp.fill(0);
            gap_add_to_bitset(gap.written(), len, tmp)?;
            return self.store.or_bits(nb, tmp, true);
        }
        if !exists {
            if len == 0 {
                self.store.put(nb, Slot::Poison(Rule::GapLength));
                return Ok(());
            }
            let mut block = vec![0; len + 1];
            // `set_gap_level(&gap_head, level)`: the level of `len + 1` entries (not past 1276)
            block[0] = with_level(head, gap_calc_level(len + 1).unwrap_or(3));
            if len > 1 {
                self.bytes.u16s(&mut block[1..], (len - 1) as u32)?;
            }
            block[len] = 65535;
            self.store.put(nb, Slot::gap(block));
            return Ok(());
        }
        if len == 0 {
            return Err(refuse(Rule::GapLength));
        }
        gap.v[0] = head;
        self.bytes.u16s(&mut gap.v[1..], (len - 1) as u32)?;
        gap.v[len] = 65535;
        gap.wrote(len + 1);
        self.store.or_gap(nb, self.s.gap.written(), true)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:4844-4934 (the gap_bienc_v3
    // case of `read_gap_block`), :4701-4739 (`decode_mins`, `decode_min_max`)
    // ```c++
    //             unsigned char head_v3 = (unsigned char) bin.get_bits(8);
    //             gap_head = (bm::gap_word_t) bin.delta16s();
    //             unsigned len = (gap_head >> 3);
    //             dst_block[0] = gap_head & bm::gap_word_t(~(3 << 1)); // clear the flags
    //             unsigned ex0_cnt{0};
    //             if (head_v3 & bm::h3f_exceptions)
    //             {
    //                 if (len < bm::gap_len_cut_off_v3)
    //                 {
    //                     if (len > 1)
    //                     {
    //                         dst_block[1] = (bm::gap_word_t) bin.delta16s();
    //                         for (unsigned k = 2; k < len; ++k)
    //                         {
    //                             auto delta = bin.delta16s();
    //                             dst_block[k] = bm::gap_word_t(dst_block[k-1] + delta);
    //                         } // for k
    //                     }
    //                     dst_block[len] = bm::gap_max_bits - 1;
    //                 }
    //                 else
    //                 {
    //                     bm::decode_min_max(bin, gap_head, head_v3, min_v, max_v);
    //                     dst_block[1] = min_v;
    //                     bm::decode_mins(bin, head_v3, min0, min1);
    //                     bin.bic_decode_u16(&dst_block[2], len-3, min_v+1, max_v);
    //                     dst_block[len-1] = max_v+1;
    //                     dst_block[len] = bm::gap_max_bits - 1;
    //                     bm::gap_restore_mins(dst_block, min0, min1);
    //                 }
    //                 unsigned h3f =
    //                     bin.decode_array(&ex0_arr[0], this->tb_wflags_, &ex0_cnt);
    //                 bool restore_v = (h3f & bm::h3f_ex_arr_1);
    //                 for (unsigned k = 0; k < ex0_cnt; ++k)
    //                     bm::gap_set_value(restore_v, dst_block, ex0_arr[k]);
    //                 if (h3f & bm::h3f_ex_arr_ex_EOC) // end-of-chain
    //                 {}
    //                 else
    //                 {
    //                     h3f =
    //                         bin.decode_array(&ex0_arr[0], this->tb_wflags_, &ex0_cnt);
    //                     restore_v = (h3f & bm::h3f_ex_arr_1);
    //                     for (unsigned k = 0; k < ex0_cnt; ++k)
    //                         bm::gap_set_value(restore_v, dst_block, ex0_arr[k]);
    //                 }
    //             } // exceptions
    //             else
    //             {
    //                 // no-exeptions decode
    //                 {
    //                     bm::decode_min_max(bin, gap_head, head_v3, min_v, max_v);
    //                     dst_block[1] = min_v;
    //                     bm::decode_mins(bin, head_v3, min0,min1);
    //                     bin.bic_decode_u16(&dst_block[2], len-3, min_v+1, max_v);
    //                     dst_block[len-1] = max_v+1;
    //                     dst_block[len] = bm::gap_max_bits - 1;
    //                 }
    //                 if (min0 || min1)
    //                     bm::gap_restore_mins(dst_block, min0, min1);
    //             }
    // ```
    /// Token 61. With exceptions, `gap_set_value` sets or clears the GAP block in
    /// `gap_temp_block_` bit by bit, as NCBI's three builds do (`gap_set_value`: on a GAP block
    /// that is not one they may differ), then the block is ORed in as the other GAP tokens.
    fn gap_v3(&mut self, nb: u32) -> Result<(), Fail> {
        let s = &mut self.s;
        let gap = &mut s.gap;
        let mut bin = BitIn::new(&mut self.bytes);
        let h3 = bin.get_bits(8)?;
        let head = bin.delta16s()? as u16;
        let len = usize::from(head >> 3);
        gap.v[0] = head & !6;
        let exceptions = h3 & 0x80 != 0;
        if exceptions && len < 4 {
            if len > 1 {
                gap.v[1] = bin.delta16s()? as u16;
                for k in 2..len {
                    gap.v[k] = u32::from(gap.v[k - 1]).wrapping_add(bin.delta16s()?) as u16;
                }
            }
            // with `len` 0 this is the header
            gap.v[len] = 65535;
            gap.wrote(len + 1);
        } else {
            // decode_min_max
            let min = if head & 2 != 0 {
                bin.get_bits(8)?
            } else {
                bin.get_16_no()?
            } as u16;
            let tail = if head & 4 != 0 {
                ((bin.get_bits(8)? as u16) << 3) | (h3 & 7) as u16
            } else {
                bin.get_16_no()? as u16
            };
            let max = 65535 - tail;
            gap.v[1] = min;
            // decode_mins
            let mut mins = [0u16; 2];
            for (min, (skip, eight)) in mins.iter_mut().zip([(0x08, 0x10), (0x40, 0x20)]) {
                if h3 & skip == 0 {
                    *min = if h3 & eight != 0 {
                        bin.gamma8()?
                    } else {
                        bin.delta16()?
                    } as u16;
                }
            }
            bin.bic_u16(
                &mut gap.v[2..],
                (len as u32).wrapping_sub(3),
                min.wrapping_add(1),
                max,
            )?;
            gap.v[len - 1] = max.wrapping_add(1);
            gap.v[len] = 65535;
            gap.wrote(len + 1);
            if exceptions || mins != [0, 0] {
                gap_restore_mins(&mut gap.v, mins[0], mins[1]);
            }
        }
        if exceptions {
            for _ in 0..2 {
                let (h3f, count) = bin.decode_array(&mut s.ex0, &mut s.wflags, 0)?;
                let value = h3f & 0x10 != 0;
                for &pos in &s.ex0.v[..count as usize] {
                    gap_set_value(gap, value, pos)?;
                }
                // `h3f_ex_arr_ex_EOC`: no second list
                if h3f & 0x20 != 0 {
                    break;
                }
            }
        }
        self.store.or_gap(nb, self.s.gap.written(), true)
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5319-5386
    // (`decode_bit_block`), :5470-5504 (`decode_arrbit`), :4280-4418 (`read_bic_arr`),
    // :4608-4629 (`read_bic_gap`), :4631-4669 (`read_digest0_block`); quoted at each form below.
    // The forms decode into the block or into `temp_block_`, ORed into the block (the inverted
    // forms invert it first). A null block becomes an empty bitmap and a GAP block a bitmap
    // first; a full block stays, and the forms that would write into it read their data without
    // a write (`IS_VALID_ADDR(blk)` false), except the inverted ones, which make it a bitmap of
    // ones (`deoptimize_block`) and decode into `temp_block_`.
    // ```c++
    //     if (!blk)
    //     {
    //         blk = bman.get_allocator().alloc_bit_block();
    //         bman.set_block(nb, blk);
    //         bm::bit_block_set(blk, 0);
    //     }
    //     else
    //     if (BM_IS_GAP(blk))
    //         blk = bman.deoptimize_block(nb);
    //     ...
    //     case bm::set_block_arr_bienc_inv:
    //     case bm::set_block_arr_bienc_inv_v3:
    //     case bm::set_block_arr_bienc_inv_v3s:
    //         if (IS_FULL_BLOCK(blk))
    //             blk = bman.deoptimize_block(nb);
    //         bm::bit_block_set(temp_block_, 0);
    //         this->read_bic_arr(dec, temp_block_, btype);
    //         bm::bit_invert(temp_block_);
    //         bm::bit_block_or(blk, temp_block_);
    //         break;
    // ```
    fn bit_block(&mut self, token: u8, nb: u32) -> Result<(), Fail> {
        let full = self.store.bit_target(nb)?;
        let inverted = matches!(token, 30 | 32 | 64 | 66);
        let dry = full && !inverted;
        let s = &mut self.s;
        let tmp = &mut s.tmp;
        tmp.fill(0);
        match token {
            30 => {
                // ```c++
                //     bm::gap_word_t len = dec.get_16();
                //     ...
                //     for (unsigned k = 0; k < len; ++k)
                //     {
                //         gap_word_t bit_idx = dec.get_16();
                //         bm::set_bit(blk, bit_idx);
                //     }
                // ```
                // (30: `bit_block_set(temp_block_, ~0u)`, `clear_bit` for each, ORed.)
                invert(tmp);
                let n = self.bytes.u16()?;
                for _ in 0..n {
                    clear_bit(tmp, u32::from(self.bytes.u16()?));
                }
            }
            31 | 32 | 57 => {
                // ```c++
                //     case bm::set_block_arr_bienc:
                //     case bm::set_block_arr_bienc_inv:
                //         min_v = dec.get_16();
                //         max_v = dec.get_16();
                //         break;
                //     case bm::set_block_arr_bienc_8bh:
                //         min_v = dec.get_8();
                //         max_delta = dec.get_8();
                //         max_v = bm::gap_word_t(65536 - max_delta);
                //         break;
                //     ...
                //     arr_len = dec.get_16();
                //     bit_in_type bin(dec);
                //     if (!IS_VALID_ADDR(blk))
                //     {
                //         bin.bic_decode_u16_dry(arr_len-2, min_v, max_v);
                //         return;
                //     }
                //     bm::set_bit(blk, min_v);
                //     bm::set_bit(blk, max_v);
                //     bin.bic_decode_u16_bitset(blk, arr_len-2, min_v, max_v);
                // ```
                let (min, max) = if token == 57 {
                    let min = u16::from(self.bytes.u8()?);
                    let delta = u32::from(self.bytes.u8()?);
                    (min, (65536 - delta) as u16)
                } else {
                    (self.bytes.u16()?, self.bytes.u16()?)
                };
                let arr_len = self.bytes.u16()?;
                if !dry {
                    set_bit(tmp, u32::from(min));
                    set_bit(tmp, u32::from(max));
                }
                BitIn::new(&mut self.bytes).bic_u16_bitset(
                    tmp,
                    u32::from(arr_len).wrapping_sub(2),
                    min,
                    max,
                    dry,
                )?;
            }
            65 | 66 => {
                // ```c++
                //             arr_len = (bm::gap_word_t) bin.delta16s();
                //             bool need_min_max = bin.get_bits(1);
                //             if (need_min_max)
                //             {
                //                 min_v = (bm::gap_word_t) bin.delta16s();
                //                 arr_len -= 2;
                //                 max_delta = (bm::gap_word_t) bin.delta16s();
                //                 max_v = bm::gap_word_t(65536 - max_delta);
                //                 if (IS_VALID_ADDR(blk))
                //                 {
                //                     bm::set_bit(tblk, min_v); bm::set_bit(tblk, max_v);
                //                 }
                //                 min_v++; max_v--;
                //             }
                //             else
                //             {
                //                 min_v = 0; max_v = 65535;
                //             }
                //             if (IS_VALID_ADDR(blk))
                //             {
                //                 bin.bic_decode_u16_bitset(tblk, arr_len, min_v, max_v);
                //             }
                //             else
                //             {
                //                 bin.bic_decode_u16_dry(arr_len, min_v, max_v);
                //             }
                // ```
                let mut bin = BitIn::new(&mut self.bytes);
                let mut arr_len = bin.delta16s()? as u16;
                let (mut min, mut max) = (0u16, 65535u16);
                if bin.get_bits(1)? != 0 {
                    min = bin.delta16s()? as u16;
                    arr_len = arr_len.wrapping_sub(2);
                    let delta = u32::from(bin.delta16s()? as u16);
                    max = (65536 - delta) as u16;
                    if !dry {
                        set_bit(tmp, u32::from(min));
                        set_bit(tmp, u32::from(max));
                    }
                    min = min.wrapping_add(1);
                    max = max.wrapping_sub(1);
                }
                bin.bic_u16_bitset(tmp, u32::from(arr_len), min, max, dry)?;
            }
            63 | 64 => {
                // ```c++
                //             unsigned h3f =
                //                 bin.decode_array(arr_s, this->tb_wflags_, &s_cnt);
                //             if (IS_VALID_ADDR(blk))
                //                 bm::bit_block_rle_set(blk, arr_s,
                //                         (bm::gap_word_t*)0, (bm::gap_word_t*)0, s_cnt, 0u);
                //             if (bool no_gap_EOC = !(h3f & bm::h3f_ex_arr_ex_EOC); no_gap_EOC)
                //             {
                //                 ...
                //                 h3f =
                //                     bin.decode_array(arr_r, this->tb_wflags_, &r_cnt);
                //                 h3f =
                //                     bin.decode_array(arr_rl, this->tb_wflags_, &s_cnt, r_cnt);
                //                 if (IS_VALID_ADDR(blk))
                //                 {
                //                     if ((h3f & bm::h3f_ex_upper2) == 1u) // BIC coder used
                //                     {
                //                         for (unsigned i = 0; i < r_cnt; ++i)
                //                             arr_rl[i] -= arr_r[i];
                //                     }
                //                     bm::bit_block_rle_set(
                //                         blk, (bm::gap_word_t*)0, arr_r, arr_rl, 0, r_cnt);
                //                 } // if VALID(blk)
                //             }
                // ```
                // NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:10069-10084
                // (`bit_block_rle_set`: `or_bit_block(blk, r[i], rlen[i]+1)`, past the block when
                // the run passes bit 65535)
                // The run lengths are `arr_rl[0..r_cnt]`, which the third array may not have
                // written (its own count is 0 or 1 for a single value): entries of an earlier
                // token, or memory NCBI did not write. NCBI writes the lengths back into
                // `ex1_arr_`.
                let mut bin = BitIn::new(&mut self.bytes);
                let (h3, s_cnt) = bin.decode_array(&mut s.ex0, &mut s.wflags, 0)?;
                if !dry {
                    for &bit in &s.ex0.v[..s_cnt as usize] {
                        set_bit(tmp, u32::from(bit));
                    }
                }
                if h3 & 0x20 == 0 {
                    let (_, r_cnt) = bin.decode_array(&mut s.ex0, &mut s.wflags, 0)?;
                    let (h3, _) = bin.decode_array(&mut s.ex1, &mut s.wflags, r_cnt)?;
                    let r_cnt = r_cnt as usize;
                    if !dry {
                        s.ex1.check(r_cnt)?;
                        let (starts, lengths) = (&s.ex0.v[..r_cnt], &mut s.ex1.v[..r_cnt]);
                        for (&start, length) in starts.iter().zip(lengths.iter_mut()) {
                            if h3 & 3 == 1 {
                                *length = length.wrapping_sub(start);
                            }
                            let end = u32::from(start) + u32::from(*length);
                            if end >= BLOCK_BITS {
                                return Err(refuse(Rule::RunRange));
                            }
                            fill_range(tmp, u32::from(start), end);
                        }
                    }
                }
            }
            33 => {
                // ```c++
                //     bm::gap_word_t head = dec.get_8();
                //     unsigned arr_len = dec.get_16();
                //     bm::gap_word_t min_v = dec.get_16();
                //     id_array_[0] = head;
                //     id_array_[1] = min_v;
                //     id_array_[arr_len] = 65535;
                //     bit_in_type bin(dec);
                //     bin.bic_decode_u16(&id_array_[2], arr_len-2, min_v, 65535);
                //     if (!IS_VALID_ADDR(blk))
                //         return;
                //     bm::gap_add_to_bitset(blk, id_array_, arr_len);
                // ```
                let head = u16::from(self.bytes.u8()?);
                let arr_len = usize::from(self.bytes.u16()?);
                let min = self.bytes.u16()?;
                let ids = &mut s.ids;
                ids.v[0] = head;
                ids.v[1] = min;
                ids.v[arr_len] = 65535;
                BitIn::new(&mut self.bytes).bic_u16(
                    &mut ids.v[2..],
                    (arr_len as u32).wrapping_sub(2),
                    min,
                    65535,
                )?;
                ids.wrote(arr_len + 1);
                if !dry {
                    gap_add_to_bitset(ids.written(), arr_len, tmp)?;
                }
            }
            _ => {
                // 34
                // ```c++
                //     bm::id64_t d0 = dec.get_64();
                //     while (d0)
                //     {
                //         bm::id64_t t = bm::bmi_blsi_u64(d0); // d & -d;
                //         unsigned wave = bm::word_bitcount64(t - 1);
                //         unsigned off = wave * bm::set_block_digest_wave_size;
                //         ...
                //                 block[off+j+0] |= dec.get_32();
                //                 block[off+j+1] |= dec.get_32();
                //                 block[off+j+2] |= dec.get_32();
                //                 block[off+j+3] |= dec.get_32();
                //         ...
                //         d0 = bm::bmi_bslr_u64(d0); // d &= d - 1;
                //     } // while
                // ```
                let mut d0 = self.bytes.u64()?;
                while d0 != 0 {
                    let wave = d0.trailing_zeros() as usize;
                    self.bytes.words(tmp, wave * 32, 32)?;
                    d0 &= d0 - 1;
                }
            }
        }
        if matches!(token, 32 | 64 | 66) {
            invert(tmp);
        }
        if !dry {
            self.store.or_into_bits(nb, tmp);
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5470-5504 (`decode_arrbit`)
    // ```c++
    //     bm::gap_word_t len = dec.get_16();
    //     if (BM_IS_GAP(blk))
    //     {
    //         blk = bman.deoptimize_block(nb);
    //     }
    //     else
    //     {
    //         if (!blk)  // block does not exists yet
    //         {
    //             ...
    //         }
    //         else
    //             if (IS_FULL_BLOCK(blk)) // nothing to do
    //             {
    //                 for (unsigned k = 0; k < len; ++k) // dry read
    //                     dec.get_16();
    //                 return;
    //             }
    //     }
    //     for (unsigned k = 0; k < len; ++k)
    //     {
    //         gap_word_t bit_idx = dec.get_16();
    //         bm::set_bit(blk, bit_idx);
    //     }
    // ```
    /// Token 16.
    fn arr_bit(&mut self, nb: u32) -> Result<(), Fail> {
        let n = self.bytes.u16()?;
        let full = self.store.bit_target(nb)?;
        let tmp = &mut self.s.tmp;
        tmp.fill(0);
        for _ in 0..n {
            set_bit(tmp, u32::from(self.bytes.u16()?));
        }
        if !full {
            self.store.or_into_bits(nb, tmp);
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5388-5420
    // (`decode_arr_sblock`: `set_bit_no_check(sb * 2^24 + arr[i])`, which leaves the size),
    // :4420-4593 (`read_bic_sb_arr`)
    // ```c++
    //     case bm::set_sblock_bienc:
    //         {
    //             sb_flag = dec.get_8();
    //             if (sb_flag & bm::sblock_flag_sb32)
    //                 *sb_idx = dec.get_32();
    //             else if (sb_flag & bm::sblock_flag_sb16)
    //                 *sb_idx = dec.get_16();
    //             else
    //                 *sb_idx = dec.get_8();
    //             if (sb_flag & bm::sblock_flag_len16)
    //                 len = dec.get_16();
    //             else
    //                 len = dec.get_8();
    //             if (!len) // there is a known issue in older version which could produce this
    //             {
    //                     throw std::logic_error(err_msg());
    //             }
    //             ... (min_v, max_v by the flags; max_v = bm::set_sub_total_bits - max_v)
    //             unsigned min0 = 0;
    //             if (sb_flag & bm::sblock_flag_dr_min)
    //             {
    //                 bool use_gamma = bin.get_bit();
    //                 if (use_gamma)
    //                     min0 = bin.gamma();
    //                 else
    //                     min0 = bin.get_16_no();
    //             }
    //             dst_arr[0] = min_v;
    //             dst_arr[len-1] = max_v;
    //             bin.bic_decode_u32_cm(&dst_arr[1], len-2, min_v, max_v);
    //             if (min0)
    //             {
    //                 bm::arr_restore_min(&dst_arr[0], len, min0);
    //             }
    //         }
    //         break;
    //     case bm::set_sblock_bienc_v3:
    //         {
    //             sb_flag = bin.get_bits(8);
    //             len = (sb_flag & bm::sblock_flag_len16) ?
    //                             bin.delta16():
    //                             bin.get_bits(8);
    //             ... (min_v, max_v, min0 by the flags, then the super-block's index)
    //             dst_arr[0] = min_v;
    //             dst_arr[len-1] = max_v;
    //             bin.bic_decode_u32_cm(&dst_arr[1], len-2, min_v+1, max_v-1);
    //             if (min0)
    //             {
    //                 bm::arr_restore_min(&dst_arr[0], len, min0);
    //             }
    //         }
    // ```
    fn super_block(&mut self, token: u8) -> Result<(), Fail> {
        const SUB_TOTAL_BITS: u32 = 1 << 24;
        let arr = &mut self.s.sb[..];
        let mut bin = BitIn::new(&mut self.bytes);
        let (len, sb, min, max, min0);
        if token == 56 {
            let flag = bin.src.u8()?;
            sb = if flag & 2 != 0 {
                bin.src.u32()?
            } else if flag & 1 != 0 {
                u32::from(bin.src.u16()?)
            } else {
                u32::from(bin.src.u8()?)
            };
            len = if flag & 0x10 != 0 {
                u32::from(bin.src.u16()?)
            } else {
                u32::from(bin.src.u8()?)
            };
            if len == 0 {
                return Err(Fail::Invalid);
            }
            let mut sized = |b24: u8, b16: u8| -> Result<u32, Fail> {
                Ok(match (flag & b24 != 0, flag & b16 != 0) {
                    (true, true) => bin.src.u32()?,
                    (true, false) => bin.src.u24()?,
                    (false, true) => u32::from(bin.src.u16()?),
                    (false, false) => u32::from(bin.src.u8()?),
                })
            };
            min = sized(0x08, 0x04)?;
            max = SUB_TOTAL_BITS.wrapping_sub(sized(0x40, 0x20)?);
            min0 = if flag & 0x80 == 0 {
                0
            } else if bin.get_bit()? != 0 {
                bin.gamma()?
            } else {
                bin.get_16_no()?
            };
        } else {
            let flag = bin.get_bits(8)?;
            len = if flag & 0x10 != 0 {
                bin.delta16()?
            } else {
                bin.get_bits(8)?
            };
            min = match (flag & 0x08 != 0, flag & 0x04 != 0) {
                (true, true) => 0,
                (true, false) => {
                    let j = bin.gamma()?;
                    j.wrapping_mul(65536).wrapping_add(bin.get_16_no()?)
                }
                (false, true) => bin.get_16_no()?,
                (false, false) => bin.get_bits(8)?,
            };
            let raw_max = match (flag & 0x40 != 0, flag & 0x20 != 0) {
                (true, true) => 0,
                (true, false) => bin.get_24_no()?,
                (false, true) => bin.get_16_no()?,
                (false, false) => bin.get_bits(8)?,
            };
            max = SUB_TOTAL_BITS.wrapping_sub(raw_max);
            min0 = if flag & 0x80 == 0 {
                0
            } else {
                match bin.gamma()? {
                    1 => bin.gamma()?,
                    2 => bin.get_bits(8)?,
                    3 => bin.get_16_no()?,
                    _ => 0,
                }
            };
            sb = if flag & 3 == 3 {
                bin.gamma()?.wrapping_sub(1)
            } else if flag & 2 != 0 {
                bin.get_32_no()?
            } else if flag & 1 != 0 {
                bin.get_16_no()?
            } else {
                bin.get_bits(8)?
            };
        }
        // `dst_arr[len-1]` with `len` 0 (68) writes 2^32 - 1 entries away.
        let last = (len as usize)
            .checked_sub(1)
            .filter(|&last| last < arr.len())
            .ok_or(refuse(Rule::Bounds))?;
        arr[0] = min;
        arr[last] = max;
        let (lo, hi) = if token == 56 {
            (min, max)
        } else {
            (min.wrapping_add(1), max.wrapping_sub(1))
        };
        bin.bic_u32(&mut arr[1..], len.wrapping_sub(2), lo, hi)?;
        let arr = &mut arr[..len as usize];
        if min0 != 0 {
            arr_restore_min_u32(arr, min0);
        }
        let from = sb.wrapping_mul(SUB_TOTAL_BITS);
        for &v in arr.iter() {
            let bit = from.wrapping_add(v);
            self.store.set_bit(bit >> 16, bit as u16)?;
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:4200-4278 (`read_id_list`)
// ```c++
//     case set_block_arrgap:
//     case set_block_arrgap_inv:
//         len = decoder.get_16();
//         decoder.get_16(dst_arr, len);
// 		break;
//     case set_block_arrgap_egamma:
//     case set_block_arrgap_egamma_inv:
//         {
//             bit_in_type bin(decoder);
//             len = (gap_word_t)bin.gamma();
//             gap_word_t prev = 0;
//             for (gap_word_t k = 0; k < len; ++k)
//             {
//                 gap_word_t bit_idx = (gap_word_t)bin.gamma();
//                 if (k == 0) --bit_idx; // TODO: optimization
//                 bit_idx = (gap_word_t)(bit_idx + prev);
//                 prev = bit_idx;
// 				dst_arr[k] = bit_idx;
//             } // for
//         }
//         break;
//     case set_block_arrgap_bienc:
//     case set_block_arrgap_bienc_inv:
//         {
//         min_v = decoder.get_16();
//         max_v = decoder.get_16();
//         bit_in_type bin(decoder);
//         len = bm::gap_word_t(bin.gamma() + 4);
//         dst_arr[0] = min_v;
//         dst_arr[len-1] = max_v;
//         bin.bic_decode_u16(&dst_arr[1], len-2, min_v, max_v);
//         break;
//         }
//     case set_block_arrgap_bienc_v2:
//     case set_block_arrgap_bienc_inv_v2:
//         {
//         len = decoder.get_16();
//         min_v = (len & 1) ?
//               decoder.get_8()
//             : decoder.get_16();
//         max_v = (len & (1<<1)) ?
//               decoder.get_8()
//             : decoder.get_16();
//         max_v = bm::gap_word_t(min_v + max_v);
//         len = bm::gap_word_t(len >> 2);
//         bit_in_type bin(decoder);
//         dst_arr[0] = min_v;
//         dst_arr[len-1] = max_v;
//         bin.bic_decode_u16(&dst_arr[1], len-2, min_v, max_v);
//         break;
//         }
// ```
/// `read_id_list(decoder, token, ids)`: the ids' count. The big-endian `get_16` of token 18
/// writes `ids[0]` also for a count of 0.
fn read_id_list(bytes: &mut Bytes<'_>, token: u8, ids: &mut Arr) -> Result<usize, Fail> {
    match token {
        18 => {
            let len = usize::from(bytes.u16()?);
            bytes.u16s(&mut ids.v, len as u32)?;
            ids.wrote(if bytes.beq { len.max(1) } else { len });
            Ok(len)
        }
        21 | 23 => {
            let mut bin = BitIn::new(bytes);
            let len = usize::from(bin.gamma()? as u16);
            let mut prev: u16 = 0;
            for (k, id) in ids.v.iter_mut().take(len).enumerate() {
                let mut bit = bin.gamma()? as u16;
                if k == 0 {
                    bit = bit.wrapping_sub(1);
                }
                prev = bit.wrapping_add(prev);
                *id = prev;
            }
            ids.wrote(len);
            Ok(len)
        }
        28 | 29 => {
            let min = bytes.u16()?;
            let max = bytes.u16()?;
            let mut bin = BitIn::new(bytes);
            let len = bin.gamma()?.wrapping_add(4) as u16;
            bic_id_list(&mut bin, ids, len, min, max)
        }
        _ => {
            // 44, 45
            let l = bytes.u16()?;
            let min = if l & 1 != 0 {
                u16::from(bytes.u8()?)
            } else {
                bytes.u16()?
            };
            let delta = if l & 2 != 0 {
                u16::from(bytes.u8()?)
            } else {
                bytes.u16()?
            };
            let mut bin = BitIn::new(bytes);
            bic_id_list(&mut bin, ids, l >> 2, min, min.wrapping_add(delta))
        }
    }
}

/// `dst_arr[0] = min_v; dst_arr[len-1] = max_v; bin.bic_decode_u16(&dst_arr[1], len-2, min_v,
/// max_v);` of `read_id_list` (a length of 0 writes before the array).
fn bic_id_list(
    bin: &mut BitIn<'_, '_>,
    ids: &mut Arr,
    len: u16,
    min: u16,
    max: u16,
) -> Result<usize, Fail> {
    let n = usize::from(len);
    let last = n.checked_sub(1).ok_or(refuse(Rule::Bounds))?;
    ids.v[0] = min;
    ids.v[last] = max;
    bin.bic_u16(&mut ids.v[1..], u32::from(len).wrapping_sub(2), min, max)?;
    ids.wrote(n);
    Ok(n)
}

#[cfg(test)]
#[path = "bitmagic_tests.rs"]
mod tests;
