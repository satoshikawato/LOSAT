//! BitMagic 8.0.1's deserialization as NCBI BLAST+ 2.17.0's object manager runs it on the data of
//! a Seq-table's `BVector-data` (line RG Q-BM, `bm/SPEC.md`): `bm::deserialize` into a
//! `bm::bvector<>` made with the `size` member, in NCBI's 32-bit BitMagic build (`BM64ADDR` not
//! defined; `BM_ASSERT` compiled out). NCBI passes the data without its length; LOSAT decodes the
//! same tokens in NCBI's order and stops at the first event of NCBI's that it does not port
//! (`Fail::Refuse`): a read past the data, an XOR reference token, a step that NCBI runs with
//! undefined behaviour or a write past its buffers, or a block whose result would depend on how
//! BitMagic holds it (runs out of order, a GAP header that wraps, a revisited block). NCBI's own
//! checks throw `std::logic_error("BM::Invalid serialization format")` (`Fail::Invalid`).
//!
//! Every token ORs bits into the vector (`bm/SPEC.md` §3.4), so the decoder keeps the set of bits
//! and the size, which is all that NCBI's index asks of the vector (`BitVector`: the size, a bit,
//! and the count of set bits up to a bit; §1.3), not BitMagic's representation of its blocks.
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
/// The largest length a GAP header's 13 bits hold (`gap_length`, `*buf >> 3`).
const GAP_LEN_MAX: u32 = 8191;
/// A block of more runs than this is kept as a bitmap.
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
    /// A GAP block of length 0, or whose length leaves the header's 13 bits (U7, U9).
    GapLength,
    /// A GAP block whose ends are not increasing below 65535 (U12).
    GapOrder,
    /// An id list that is empty or not increasing (U9).
    IdList,
    /// A window of 0 words (a division by zero) or wider than its array (U11).
    Window,
    /// Fewer run lengths than runs (NCBI reads the rest from memory that it did not write).
    RunLengths,
    /// A run of bits past the block (U10).
    RunRange,
    /// A big-endian GAP block of length 1 in a block that may already exist: NCBI then reads two
    /// bytes more or not, by how BitMagic holds the block.
    Revisit,
}

/// The `what()` of NCBI's `std::logic_error` (`Fail::Invalid`).
pub(crate) const INVALID: &str = "BM::Invalid serialization format";

const fn refuse(rule: Rule) -> Fail {
    Fail::Refuse(rule)
}

/// The data's size and set bits as NCBI's index asks them (`bm::bvector<>::size`, `get_bit`,
/// `count_range`), kept by blocks: runs of full blocks, blocks as runs of set bits, blocks as
/// bitmaps.
#[derive(Debug)]
pub(crate) struct BitVector {
    size: u32,
    segs: Vec<Seg>,
    count: u64,
}

#[derive(Debug)]
struct Seg {
    first: u32,
    last: u32,
    /// The set bits of the segments before.
    before: u64,
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

    /// `bvector<>::get_bit(n)` (bm.h:3604-3631).
    pub(crate) fn get_bit(&self, n: u32) -> bool {
        let p = n & (BLOCK_BITS - 1);
        match self.seg(n >> 16).1.map(|seg| &seg.kind) {
            None => false,
            Some(Kind::Full) => true,
            Some(Kind::Runs(runs, _)) => {
                let at = runs.partition_point(|run| u32::from(run.1) < p);
                runs.get(at).is_some_and(|run| u32::from(run.0) <= p)
            }
            Some(Kind::Bits(words, _)) => words[(p >> 6) as usize & 1023] >> (p & 63) & 1 != 0,
        }
    }

    /// `bvector<>::count_range(0, n)` (bm.h:3244-3254): the set bits of `[0, n]`.
    pub(crate) fn count_to(&self, n: u32) -> u64 {
        let block = n >> 16;
        let p = n & (BLOCK_BITS - 1);
        let (at, seg) = self.seg(block);
        let Some(seg) = seg else {
            return self.segs.get(at).map_or(self.count, |seg| seg.before);
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
        };
        seg.before + inside
    }

    /// All set bits.
    #[cfg(test)]
    pub(crate) fn count(&self) -> u64 {
        self.count
    }

    /// The runs of set bits `(first, last)` in order.
    #[cfg(test)]
    pub(crate) fn ranges(&self) -> Vec<(u32, u32)> {
        let mut out: Vec<(u32, u32)> = Vec::new();
        let mut push = |a: u32, b: u32| match out.last_mut() {
            Some(last) if last.1.checked_add(1) == Some(a) => last.1 = b,
            _ => out.push((a, b)),
        };
        for seg in &self.segs {
            match &seg.kind {
                Kind::Full => push(seg.first << 16, (seg.last << 16) | 0xFFFF),
                Kind::Runs(runs, _) => {
                    for &(a, b) in runs.iter() {
                        push(
                            seg.first << 16 | u32::from(a),
                            seg.first << 16 | u32::from(b),
                        );
                    }
                }
                Kind::Bits(words, _) => {
                    for (a, b) in block_runs(words, usize::MAX).unwrap_or_default() {
                        push(
                            seg.first << 16 | u32::from(a),
                            seg.first << 16 | u32::from(b),
                        );
                    }
                }
            }
        }
        out
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
    Ok((decoder.store.finish(decoder.size), end))
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
    // :1219-1236 (big-endian; `count` is at least 1 at every caller here)
    /// `get_32(block + first, count)` into a cleared block.
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
    /// block's 65536 bits writes past it, U4). NCBI's dry variant for a full block reads the same
    /// bits; LOSAT checks its values as this one's (LOSAT does not know whether NCBI's block is
    /// full; a full block keeps its bits either way).
    fn bic_u16_bitset(&mut self, block: &mut Block, sz: u32, lo: u16, hi: u16) -> Result<(), Fail> {
        if sz == 0 {
            return Ok(());
        }
        self.bic_u16_cm_bitset(block, sz, lo, hi)
    }

    fn bic_u16_cm_bitset(
        &mut self,
        block: &mut Block,
        mut sz: u32,
        mut lo: u16,
        hi: u16,
    ) -> Result<(), Fail> {
        loop {
            let r = u32::from(hi)
                .wrapping_sub(u32::from(lo))
                .wrapping_sub(sz)
                .wrapping_add(1);
            let mid = sz >> 1;
            let val = self
                .bic_value(r)?
                .wrapping_add(u32::from(lo))
                .wrapping_add(mid);
            if val >= BLOCK_BITS {
                return Err(refuse(Rule::BicWrite));
            }
            set_bit(block, val);
            if sz <= 1 {
                return Ok(());
            }
            self.bic_u16_cm_bitset(block, mid, lo, val.wrapping_sub(1) as u16)?;
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
    /// `EX_ENTRIES`): the header flags and the size.
    fn decode_array(
        &mut self,
        arr: &mut [u16],
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
            *arr.first_mut().ok_or(refuse(Rule::Bounds))? = value;
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
        // Every form writes `arr[0..sz]` (and `arr[0]` with `sz` 0).
        let arr = arr
            .get_mut(..(sz as usize).max(1))
            .ok_or(refuse(Rule::Bounds))?;
        if h3 & 3 == 0 {
            // delta-gamma
            let mut prev = if h3 & 0x80 != 0 {
                0
            } else {
                self.gamma()? as u16
            };
            arr[0] = prev;
            for value in arr.iter_mut().take(sz as usize).skip(1) {
                let delta = self.gamma()?;
                prev = u32::from(prev)
                    .wrapping_add(delta)
                    .wrapping_add(u32::from(min0)) as u16;
                *value = prev;
            }
            return Ok((h3, sz));
        }
        if h3 & 2 != 0 {
            // gamma
            let zero_correct = u16::from(h3 & 0x80 != 0);
            for value in arr.iter_mut().take(sz as usize) {
                *value = (self.gamma()? as u16)
                    .wrapping_sub(zero_correct)
                    .wrapping_add(min0);
            }
            return Ok((h3, sz));
        }
        // BIC-DR
        if h3 & 0x80 != 0 {
            let mut min_v = self.get_16_no()? as u16;
            let mut max_v = self.get_16_no()? as u16;
            arr[0] = min_v;
            *arr.get_mut((sz as usize).wrapping_sub(1))
                .ok_or(refuse(Rule::Bounds))? = max_v;
            if sz == 2 {
                return Ok((h3, sz));
            }
            min_v = min_v.wrapping_add(1);
            max_v = max_v.wrapping_sub(1);
            self.bic_u16(&mut arr[1..], sz.wrapping_sub(2), min_v, max_v)?;
        } else {
            self.bic_u16(arr, sz, 0, 65535)?;
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
            self.bic_u16_bitset(wflags, wcnt, 1, max_wd as u16)?;
            arr_restore_min_w(arr, sz, win_size, min0, wflags)?;
        } else if min0 != 0 {
            arr_restore_min_u16(arr, sz, min0);
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
// first loop reads `arr[1..wlen]`, which `decode_array` wrote only up to `arr_len`: a window wider
// than the array reads memory that NCBI did not write, refused)
fn arr_restore_min_w(
    arr: &mut [u16],
    arr_len: u32,
    mut wlen: u32,
    min0: u16,
    wflags: &Block,
) -> Result<(), Fail> {
    if wlen > arr_len {
        return Err(refuse(Rule::Window));
    }
    let arr = &mut arr[..arr_len as usize];
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

fn get(block: &Block, p: u32) -> bool {
    block[(p >> 6) as usize & 1023] >> (p & 63) & 1 != 0
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

/// The set bits of a GAP block `gap[0..=len]` (`gap[0]` its header, whose lowest bit is the first
/// run's value), refused unless its ends are increasing below 65535 and it ends with 65535
/// (`gap_add_to_bitset`, `gap_operation_or`, BitMagic's searches assume it, U12).
fn gap_bits(gap: &[u16], len: usize, block: &mut Block) -> Result<(), Fail> {
    if len == 0 {
        return Err(refuse(Rule::GapLength));
    }
    let ends = gap.get(1..=len).ok_or(refuse(Rule::Bounds))?;
    let (&last, inner) = ends.split_last().ok_or(refuse(Rule::GapLength))?;
    if last != 65535 || inner.contains(&65535) || inner.windows(2).any(|pair| pair[1] <= pair[0]) {
        return Err(refuse(Rule::GapOrder));
    }
    block.fill(0);
    let mut value = gap[0] & 1 != 0;
    let mut from = 0u32;
    for &end in ends {
        if value {
            fill_range(block, from, u32::from(end));
        }
        from = u32::from(end) + 1;
        value = !value;
    }
    Ok(())
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
/// The bits that `gap_add_to_bitset(dest, gap, len)` (and `gap_convert_to_bitset`, which passes
/// the header's length) sets into a cleared block: the runs of 1s whose end is at index `len` or
/// before, with no other check on the GAP. A run that ends where it starts sets nothing
/// (`or_bit_block` with a count of 0 writes nothing, bmfunc.h:4517-4546); one that ends before it
/// starts wraps its count (a write past the block): refused.
fn gap_add_to_bitset(gap: &[u16], len: usize, block: &mut Block) -> Result<(), Fail> {
    block.fill(0);
    let mut k = 2;
    if gap[0] & 1 != 0 {
        fill_range(block, 0, u32::from(gap[1]));
        k = 3;
    }
    while k <= len {
        let (from, to) = (gap[k - 1], gap[k]);
        if to < from {
            return Err(refuse(Rule::GapOrder));
        }
        if to > from {
            fill_range(block, u32::from(from) + 1, u32::from(to));
        }
        k += 2;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/include/util/bitset/bmfunc.h:5409-5418 (`gap_calc_level` with
// the default `gap_len_table`, 128, 256, 512, 1280: -1 above 1276), called with `gap_length`, one
// more than the GAP's last index
/// Whether a GAP block of last index `len` is too long for a GAP block and converted to bits.
fn gap_converts(len: usize) -> bool {
    len + 1 > 1280 - 4
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
/// `gap_set_array(gap, ids, len)` after `gap[0] = 0` (ids in any order): the GAP's last index,
/// whose header keeps its low 13 bits.
fn gap_set_array(gap: &mut [u16], ids: &[u16]) -> Result<usize, Fail> {
    fn put(gap: &mut [u16], at: usize, value: u16) -> Result<(), Fail> {
        *gap.get_mut(at).ok_or(refuse(Rule::Bounds))? = value;
        Ok(())
    }
    let (&first, rest) = ids.split_first().ok_or(refuse(Rule::IdList))?;
    gap[0] = 1 << 3;
    let mut at = 1;
    if first != 0 {
        put(gap, at, first - 1)?;
        at += 1;
    } else {
        gap[0] += 1;
    }
    let mut prev = first;
    let mut acc = prev;
    for &curr in rest {
        if u32::from(curr) == u32::from(prev) + 1 {
            acc = acc.wrapping_add(1);
        } else {
            put(gap, at, acc)?;
            acc = curr;
            put(gap, at + 1, curr.wrapping_sub(1))?;
            at += 2;
        }
        prev = curr;
    }
    put(gap, at, acc)?;
    if acc != 65535 {
        at += 1;
        put(gap, at, 65535)?;
    }
    gap[0] = (gap[0] & 7).wrapping_add((at as u16) << 3);
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
/// `gap_add_value(gap, pos)`: sets bit `pos` of a GAP block whose last run is of 0s and holds
/// `pos` (token 24 checks its ids increasing, from an empty block, so that this holds). The GAP's
/// length, refused when it leaves the header's 13 bits.
fn gap_add_value(gap: &mut [u16], pos: u16) -> Result<usize, Fail> {
    let mut end = usize::from(gap[0] >> 3);
    let cur = end;
    let prev = end.checked_sub(1).ok_or(refuse(Rule::GapLength))?;
    if pos == 0 {
        gap[0] ^= 1;
        if gap[1] != 0 {
            gap.copy_within(1..end, 2);
            gap[1] = 0;
            end += 1;
        } else {
            end -= 1;
            gap.copy_within(2..cur.max(3), 1);
        }
    } else if cur > 1 && u32::from(gap[prev]) + 1 == u32::from(pos) {
        gap[prev] = gap[prev].wrapping_add(1);
        if gap[prev] == gap[cur] {
            end -= 1;
        }
    } else if gap[cur] == pos {
        gap[cur] = gap[cur].wrapping_sub(1);
        end += 1;
    } else {
        gap[cur] = pos.wrapping_sub(1);
        gap[cur + 1] = pos;
        end += 2;
    }
    if end > GAP_LEN_MAX as usize {
        return Err(refuse(Rule::GapLength));
    }
    gap[0] = (gap[0] & 7) + ((end as u16) << 3);
    gap[end] = 65535;
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

/// A block's bits while decoding: runs of set bits, a bitmap, or full. A block that a token wrote
/// (`Some`, possibly without bits) may exist in NCBI's vector; one that no token wrote does not.
enum Slot {
    Runs(Vec<(u16, u16)>),
    Bits(Box<Block>),
    Full,
}

impl Slot {
    fn of(block: &Block) -> Slot {
        if block.iter().all(|&word| word == u64::MAX) {
            return Slot::Full;
        }
        match block_runs(block, RUNS_MAX) {
            Some(runs) => Slot::Runs(runs),
            None => Slot::Bits(Box::new(*block)),
        }
    }
}

/// The blocks by super-block of 256 (BitMagic's two levels), and the blocks that hold a GAP block
/// that is not one (U12): NCBI stores it without looking at it, and refuses nothing until a token
/// writes into the block or the decoding ends, when LOSAT refuses the data.
struct Store {
    top: Vec<Option<Box<[Option<Slot>; 256]>>>,
    work: Box<Block>,
    group: Box<Block>,
    poisoned: Vec<u32>,
    poison: Option<Rule>,
}

type Top = [Option<Box<[Option<Slot>; 256]>>];

impl Default for Store {
    fn default() -> Self {
        Store {
            top: (0..256).map(|_| None).collect(),
            work: Box::new([0; 1024]),
            group: Box::new([0; 1024]),
            poisoned: Vec::new(),
            poison: None,
        }
    }
}

fn slot_in(top: &mut Top, nb: u32) -> &mut Option<Slot> {
    let sub = top[(nb >> 8) as usize & 255]
        .get_or_insert_with(|| Box::new(std::array::from_fn(|_| None)));
    &mut sub[nb as usize & 255]
}

/// ORs `bits` into the slot.
fn or_slot(slot: &mut Option<Slot>, work: &mut Block, bits: &Block) {
    *slot = Some(match slot.take() {
        None => Slot::of(bits),
        Some(Slot::Full) => Slot::Full,
        Some(Slot::Runs(runs)) => {
            work.copy_from_slice(bits);
            for (a, b) in runs {
                fill_range(work, u32::from(a), u32::from(b));
            }
            Slot::of(work)
        }
        Some(Slot::Bits(mut words)) => {
            for (word, bit) in words.iter_mut().zip(bits.iter()) {
                *word |= bit;
            }
            match Slot::of(&words) {
                Slot::Bits(_) => Slot::Bits(words),
                other => other,
            }
        }
    });
}

impl Store {
    /// Whether a token wrote block `nb`.
    fn touched(&self, nb: u32) -> bool {
        self.top[(nb >> 8) as usize & 255]
            .as_ref()
            .is_some_and(|sub| sub[nb as usize & 255].is_some())
    }

    fn check(&self, nb: u32) -> Result<(), Fail> {
        match self.poison {
            Some(rule) if self.poisoned.contains(&nb) => Err(refuse(rule)),
            _ => Ok(()),
        }
    }

    /// Block `nb`, which no token wrote, now holds a GAP block that is not one.
    fn poison(&mut self, nb: u32, rule: Rule) {
        *slot_in(&mut self.top, nb) = Some(Slot::Runs(Vec::new()));
        self.poisoned.push(nb);
        self.poison.get_or_insert(rule);
    }

    /// ORs `bits` into block `nb`.
    fn or_block(&mut self, nb: u32, bits: &Block) -> Result<(), Fail> {
        self.check(nb)?;
        or_slot(slot_in(&mut self.top, nb), &mut self.work, bits);
        Ok(())
    }

    fn set_full(&mut self, nb: u32) -> Result<(), Fail> {
        self.check(nb)?;
        *slot_in(&mut self.top, nb) = Some(Slot::Full);
        Ok(())
    }

    /// Sets the bits `positions` (any order).
    fn set_bits(&mut self, positions: &mut [u32]) -> Result<(), Fail> {
        positions.sort_unstable();
        for group in positions.chunk_by(|a, b| a >> 16 == b >> 16) {
            let nb = group[0] >> 16;
            self.check(nb)?;
            self.group.fill(0);
            for &p in group {
                set_bit(&mut self.group, p & (BLOCK_BITS - 1));
            }
            or_slot(slot_in(&mut self.top, nb), &mut self.work, &self.group);
        }
        Ok(())
    }

    /// Sets the bits `left..=right`.
    fn set_range(&mut self, left: u32, right: u32) -> Result<(), Fail> {
        let (first, last) = (left >> 16, right >> 16);
        for nb in first..=last {
            let from = if nb == first { left & 0xFFFF } else { 0 };
            let to = if nb == last { right & 0xFFFF } else { 0xFFFF };
            if from == 0 && to == 0xFFFF {
                self.set_full(nb)?;
            } else {
                self.check(nb)?;
                self.group.fill(0);
                fill_range(&mut self.group, from, to);
                or_slot(slot_in(&mut self.top, nb), &mut self.work, &self.group);
            }
        }
        Ok(())
    }

    fn finish(self, size: u32) -> BitVector {
        let mut segs: Vec<Seg> = Vec::new();
        let mut count = 0u64;
        for (i, sub) in self.top.into_iter().enumerate() {
            let Some(sub) = sub else { continue };
            for (j, slot) in sub.into_iter().enumerate() {
                let nb = (i as u32) << 8 | j as u32;
                let kind = match slot {
                    None => continue,
                    Some(Slot::Full) => {
                        if let Some(last) = segs.last_mut() {
                            if matches!(last.kind, Kind::Full) && last.last + 1 == nb {
                                last.last = nb;
                                count += u64::from(BLOCK_BITS);
                                continue;
                            }
                        }
                        Kind::Full
                    }
                    Some(Slot::Runs(runs)) => {
                        if runs.is_empty() {
                            continue;
                        }
                        let mut before = Vec::with_capacity(runs.len() + 1);
                        let mut sum = 0u32;
                        before.push(sum);
                        for &(a, b) in &runs {
                            sum += u32::from(b) - u32::from(a) + 1;
                            before.push(sum);
                        }
                        Kind::Runs(runs.into_boxed_slice(), before.into_boxed_slice())
                    }
                    Some(Slot::Bits(words)) => {
                        let mut chunks = Box::new([0u32; 64]);
                        let mut sum = 0u32;
                        for (chunk, words) in chunks.iter_mut().zip(words.chunks(16)) {
                            *chunk = sum;
                            sum += words.iter().map(|word| word.count_ones()).sum::<u32>();
                        }
                        Kind::Bits(words, chunks)
                    }
                };
                let bits = match &kind {
                    Kind::Full => u64::from(BLOCK_BITS),
                    Kind::Runs(_, before) => u64::from(before[before.len() - 1]),
                    Kind::Bits(words, _) => words.iter().map(|w| u64::from(w.count_ones())).sum(),
                };
                segs.push(Seg {
                    first: nb,
                    last: nb,
                    before: count,
                    kind,
                });
                count += bits;
            }
        }
        BitVector { size, segs, count }
    }
}

/// NCBI's work arrays of the deserializer (bmserial.h:5110-5129), and the block a token decodes
/// into (`temp_block_`).
struct Scratch {
    gap: Vec<u16>,
    ids: Vec<u16>,
    sb: Vec<u32>,
    ex0: Vec<u16>,
    ex1: Vec<u16>,
    wflags: Box<Block>,
    tmp: Box<Block>,
}

impl Default for Scratch {
    fn default() -> Self {
        Scratch {
            gap: vec![0; GAP_ENTRIES],
            ids: vec![0; GAP_ENTRIES],
            sb: vec![0; GAP_ENTRIES],
            ex0: vec![0; EX_ENTRIES],
            ex1: vec![0; EX_ENTRIES],
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
            self.store.set_bits(&mut ids)?;
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
                    for block in nb..TOTAL_BLOCKS {
                        self.store.set_full(block)?;
                    }
                    Go::End
                }
                2 => {
                    self.store.set_full(nb)?;
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
                    let tmp = &mut self.s.tmp;
                    tmp.fill(0);
                    self.bytes.words(tmp, 0, 2048)?;
                    self.store.or_block(nb, tmp)?;
                    Go::Next
                }
                19 => {
                    let bit = u32::from(self.bytes.u16()?).wrapping_add(nb << 16);
                    self.store.set_bits(&mut [bit])?;
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
                16 | 30 | 31 | 32 | 33 | 34 | 57 | 63 | 64 | 65 | 66 => {
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
        // a block that holds what is not a GAP block is queried by its representation
        if let Some(rule) = self.store.poison {
            return Err(refuse(rule));
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
        self.store.or_block(nb, tmp)?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5443-5468
    // (`decode_block_bit_interval`: words `head..=tail`; outside `head <= tail < 2048` NCBI
    // writes past the block or, for `tail == head - 1`, writes by how it holds the block)
    fn interval(&mut self, nb: u32) -> Result<(), Fail> {
        let head = u32::from(self.bytes.u16()?);
        let tail = u32::from(self.bytes.u16()?);
        if head > tail || tail >= 2048 {
            return Err(refuse(Rule::Interval));
        }
        let tmp = &mut self.s.tmp;
        tmp.fill(0);
        self.bytes
            .words(tmp, head as usize, (tail - head + 1) as usize)?;
        self.store.or_block(nb, tmp)?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5154-5317
    // (`deserialize_gap`), :4744-5034 (`read_gap_block`), :4200-4278 (`read_id_list`); quoted at
    // each form below. Each form makes a GAP block in `gap_temp_block_` (or, for 14 and 15 into an
    // empty block, the block itself). NCBI then converts it to bits (`gap_convert_to_bitset`: 14
    // and 15 and the id lists of `deserialize_gap` above 1276 entries, and through
    // `combine_operation_block_or` -> `clone_gap_block` a long GAP into an empty block), which
    // needs no more than runs that end after they start; otherwise it keeps or merges it as a GAP
    // block, which must be a GAP (U12). Either is ORed into block `nb`.
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
        let gap = &mut s.gap[..];
        // `len`: the GAP's length as its header holds it (`gap_length - 1`, its last index), and
        // whether NCBI converts it to bits whatever block `nb` holds.
        let (len, convert) = match token {
            14 | 15 => {
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
                //         ...
                //         if (blk == 0)
                //         {
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
                // A length of 0 writes the header's place or reads 2^32 - 1 values; a length of 1
                // in a block that exists reads one value with the big-endian decoder and none with
                // the little-endian one (`get_16(p, 0)`), so the bytes read depend on the block.
                let head = self.bytes.u16()?;
                let len = usize::from(head >> 3);
                if len == 0 {
                    // into an empty block: a GAP block whose header is 65535, written without a
                    // look; into one that exists: `get_16(p, 2^32 - 1)`
                    if self.store.touched(nb) {
                        return Err(refuse(Rule::GapLength));
                    }
                    self.store.poison(nb, Rule::GapLength);
                    return Ok(());
                }
                if len == 1 && self.bytes.beq && self.store.touched(nb) && !gap_converts(len) {
                    return Err(refuse(Rule::Revisit));
                }
                gap[0] = head;
                if len > 1 {
                    self.bytes.u16s(&mut gap[1..], (len - 1) as u32)?;
                }
                gap[len] = 65535;
                (len, gap_converts(len))
            }
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
                // ```
                // (23 and 29 in `read_gap_block`: the same, inverted, and no conversion.) The level
                // is of the GAP's length, the conversion of the header's 13 bits of it. An empty
                // list reads `arr[0]`, which no token may have written.
                let n = read_id_list(&mut self.bytes, token, &mut s.ids)?;
                let ids = &s.ids[..n];
                if ids.is_empty() {
                    return Err(refuse(Rule::IdList));
                }
                let full = gap_set_array(gap, ids)?;
                if matches!(token, 23 | 29 | 45) {
                    gap[0] ^= 1;
                }
                let convert = !matches!(token, 23 | 29) && gap_converts(full);
                (usize::from(gap[0] >> 3), convert)
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
                gap[0] = (gap[0] & 6) + (1 << 3);
                gap[1] = 65535;
                let n = self.bytes.u16()?;
                let mut prev: Option<u16> = None;
                let mut len = 1;
                for _ in 0..n {
                    let id = self.bytes.u16()?;
                    if prev.is_some_and(|prev| id <= prev) {
                        return Err(refuse(Rule::IdList));
                    }
                    prev = Some(id);
                    len = gap_add_value(gap, id)?;
                }
                gap[0] ^= 1;
                (len, false)
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
                let head = self.bytes.u16()?;
                let len = usize::from(head >> 3);
                if len == 0 {
                    return Err(refuse(Rule::GapLength));
                }
                gap[0] = head;
                let mut bin = BitIn::new(&mut self.bytes);
                let mut sum = (bin.gamma()? as u16).wrapping_sub(1);
                gap[1] = sum;
                for k in 2..len {
                    sum = sum.wrapping_add(bin.gamma()? as u16);
                    gap[k] = sum;
                }
                gap[len] = 65535;
                (len, false)
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
                gap[0] = head;
                let min = self.bytes.u16()?;
                gap[1] = min;
                BitIn::new(&mut self.bytes).bic_u16(
                    &mut gap[2..],
                    (len as u32).wrapping_sub(2),
                    min,
                    65535,
                )?;
                gap[len] = 65535;
                (len, false)
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
                gap[0] = head & !6;
                gap[1] = min;
                BitIn::new(&mut self.bytes).bic_u16(
                    &mut gap[2..],
                    (len as u32).wrapping_sub(3),
                    min,
                    max,
                )?;
                gap[len - 1] = max;
                gap[len] = 65535;
                (len, false)
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
                gap[0] = head & !6;
                gap[1] = min;
                bin.bic_u16(&mut gap[2..], (len as u32).wrapping_sub(3), min, max)?;
                gap[len - 1] = max;
                gap[len] = 65535;
                (len, false)
            }
            61 => return self.gap_v3(nb),
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
                // A length of 0 writes 65535 into the header, one above 65535 writes past
                // `gap_temp_block_`; the header keeps the low 13 bits of the length.
                let mut bin = BitIn::new(&mut self.bytes);
                let len = bin.gamma()?.wrapping_add(1);
                let start = bin.get_bit()?;
                let use_gamma = bin.get_bit()?;
                if len == 0 {
                    return Err(refuse(Rule::GapLength));
                }
                if len as usize >= GAP_ENTRIES {
                    return Err(refuse(Rule::Bounds));
                }
                let len = len as usize;
                gap[0] = ((len as u16) << 3) | start as u16;
                if use_gamma != 0 {
                    let mut prev = bin.gamma8()? as u16;
                    gap[1] = prev;
                    for k in 2..len {
                        prev = u32::from(prev).wrapping_add(bin.gamma8()?) as u16;
                        gap[k] = prev;
                    }
                } else {
                    for k in 1..len {
                        gap[k] = bin.get_16_no()? as u16;
                    }
                }
                gap[len] = 65535;
                (usize::from(gap[0] >> 3), false)
            }
        };
        let touched = self.store.touched(nb);
        if convert || (!touched && gap_converts(len)) {
            gap_add_to_bitset(&s.gap, len, &mut s.tmp)?;
        } else if let Err(fail) = gap_bits(&s.gap, len, &mut s.tmp) {
            // kept as the block, or merged with the block there
            match fail {
                Fail::Refuse(rule) if !touched => self.store.poison(nb, rule),
                _ => return Err(fail),
            }
            return Ok(());
        }
        self.store.or_block(nb, &s.tmp)?;
        Ok(())
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
    /// Token 61. The exceptions set or clear bits of the GAP block (`gap_set_value`), whose
    /// length must stay in the header's 13 bits.
    fn gap_v3(&mut self, nb: u32) -> Result<(), Fail> {
        let s = &mut self.s;
        let gap = &mut s.gap[..];
        let mut bin = BitIn::new(&mut self.bytes);
        let h3 = bin.get_bits(8)?;
        let head = bin.delta16s()? as u16;
        let len = usize::from(head >> 3);
        gap[0] = head & !6;
        let exceptions = h3 & 0x80 != 0;
        if exceptions && len < 4 {
            if len > 1 {
                gap[1] = bin.delta16s()? as u16;
                for k in 2..len {
                    gap[k] = u32::from(gap[k - 1]).wrapping_add(bin.delta16s()?) as u16;
                }
            }
            if len == 0 {
                return Err(refuse(Rule::GapLength));
            }
            gap[len] = 65535;
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
            gap[1] = min;
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
                &mut gap[2..],
                (len as u32).wrapping_sub(3),
                min.wrapping_add(1),
                max,
            )?;
            gap[len - 1] = max.wrapping_add(1);
            gap[len] = 65535;
            if exceptions || mins != [0, 0] {
                gap_restore_mins(gap, mins[0], mins[1]);
            }
        }
        let tmp = &mut s.tmp;
        // `gap_set_value` needs a GAP; without exceptions the block goes the way of the other forms.
        let touched = self.store.touched(nb);
        if !exceptions && !touched && gap_converts(len) {
            gap_add_to_bitset(gap, len, tmp)?;
        } else if let Err(fail) = gap_bits(gap, len, tmp) {
            match fail {
                Fail::Refuse(rule) if !exceptions && !touched => self.store.poison(nb, rule),
                _ => return Err(fail),
            }
            return Ok(());
        }
        if exceptions {
            // The GAP's length (`end`) as `gap_set_value` keeps it: one more than the changes of
            // value between neighbouring bits.
            let mut runs = len as u32;
            for _ in 0..2 {
                let (h3f, count) = bin.decode_array(&mut s.ex0, &mut s.wflags, 0)?;
                let value = h3f & 0x10 != 0;
                for &pos in s.ex0.iter().take(count as usize) {
                    let pos = u32::from(pos);
                    if get(tmp, pos) == value {
                        continue;
                    }
                    for q in [pos.wrapping_sub(1), pos + 1] {
                        if q < BLOCK_BITS {
                            if get(tmp, q) == value {
                                runs -= 1;
                            } else {
                                runs += 1;
                            }
                        }
                    }
                    if value {
                        set_bit(tmp, pos);
                    } else {
                        clear_bit(tmp, pos);
                    }
                    if runs > GAP_LEN_MAX {
                        return Err(refuse(Rule::GapLength));
                    }
                }
                // `h3f_ex_arr_ex_EOC`: no second list
                if h3f & 0x20 != 0 {
                    break;
                }
            }
        }
        self.store.or_block(nb, tmp)?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/util/bitset/bmserial.h:5319-5386
    // (`decode_bit_block`), :5470-5504 (`decode_arrbit`), :4280-4418 (`read_bic_arr`),
    // :4608-4629 (`read_bic_gap`), :4631-4669 (`read_digest0_block`); quoted at each form below.
    // The forms decode into the block or into `temp_block_`, ORed into the block (the inverted
    // forms invert it first); NCBI decodes into a full block without writing, which keeps its
    // bits.
    fn bit_block(&mut self, token: u8, nb: u32) -> Result<(), Fail> {
        let s = &mut self.s;
        let tmp = &mut s.tmp;
        tmp.fill(0);
        match token {
            16 | 30 => {
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
                if token == 30 {
                    invert(tmp);
                }
                let n = self.bytes.u16()?;
                for _ in 0..n {
                    let bit = u32::from(self.bytes.u16()?);
                    if token == 30 {
                        clear_bit(tmp, bit);
                    } else {
                        set_bit(tmp, bit);
                    }
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
                //     ...
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
                set_bit(tmp, u32::from(min));
                set_bit(tmp, u32::from(max));
                BitIn::new(&mut self.bytes).bic_u16_bitset(
                    tmp,
                    u32::from(arr_len).wrapping_sub(2),
                    min,
                    max,
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
                // ```
                let mut bin = BitIn::new(&mut self.bytes);
                let mut arr_len = bin.delta16s()? as u16;
                let (mut min, mut max) = (0u16, 65535u16);
                if bin.get_bits(1)? != 0 {
                    min = bin.delta16s()? as u16;
                    arr_len = arr_len.wrapping_sub(2);
                    let delta = u32::from(bin.delta16s()? as u16);
                    max = (65536 - delta) as u16;
                    set_bit(tmp, u32::from(min));
                    set_bit(tmp, u32::from(max));
                    min = min.wrapping_add(1);
                    max = max.wrapping_sub(1);
                }
                bin.bic_u16_bitset(tmp, u32::from(arr_len), min, max)?;
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
                let mut bin = BitIn::new(&mut self.bytes);
                let (h3, s_cnt) = bin.decode_array(&mut s.ex0, &mut s.wflags, 0)?;
                for &bit in s.ex0.iter().take(s_cnt as usize) {
                    set_bit(tmp, u32::from(bit));
                }
                if h3 & 0x20 == 0 {
                    let (_, r_cnt) = bin.decode_array(&mut s.ex0, &mut s.wflags, 0)?;
                    let (h3, rl_cnt) = bin.decode_array(&mut s.ex1, &mut s.wflags, r_cnt)?;
                    if rl_cnt < r_cnt {
                        return Err(refuse(Rule::RunLengths));
                    }
                    for (&start, &length) in s.ex0.iter().zip(s.ex1.iter()).take(r_cnt as usize) {
                        let length = if h3 & 3 == 1 {
                            length.wrapping_sub(start)
                        } else {
                            length
                        };
                        let end = u32::from(start) + u32::from(length);
                        if end >= BLOCK_BITS {
                            return Err(refuse(Rule::RunRange));
                        }
                        fill_range(tmp, u32::from(start), end);
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
                let ids = &mut s.ids[..];
                ids[0] = head;
                ids[1] = min;
                ids[arr_len] = 65535;
                BitIn::new(&mut self.bytes).bic_u16(
                    &mut ids[2..],
                    (arr_len as u32).wrapping_sub(2),
                    min,
                    65535,
                )?;
                // `gap_add_to_bitset(blk, id_array_, arr_len)`, not into a full block (no write)
                gap_add_to_bitset(ids, arr_len, tmp)?;
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
        self.store.or_block(nb, tmp)?;
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
        let mut bits: Vec<u32> = arr.iter().map(|&v| from.wrapping_add(v)).collect();
        self.store.set_bits(&mut bits)?;
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
/// `read_id_list(decoder, token, ids)`: the ids' count.
fn read_id_list(bytes: &mut Bytes<'_>, token: u8, ids: &mut [u16]) -> Result<usize, Fail> {
    match token {
        18 => {
            let len = bytes.u16()?;
            bytes.u16s(ids, u32::from(len))?;
            Ok(usize::from(len))
        }
        21 | 23 => {
            let mut bin = BitIn::new(bytes);
            let len = bin.gamma()? as u16;
            let mut prev: u16 = 0;
            for (k, id) in ids.iter_mut().take(usize::from(len)).enumerate() {
                let mut bit = bin.gamma()? as u16;
                if k == 0 {
                    bit = bit.wrapping_sub(1);
                }
                prev = bit.wrapping_add(prev);
                *id = prev;
            }
            Ok(usize::from(len))
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
    ids: &mut [u16],
    len: u16,
    min: u16,
    max: u16,
) -> Result<usize, Fail> {
    let n = usize::from(len);
    let last = n.checked_sub(1).ok_or(refuse(Rule::Bounds))?;
    ids[0] = min;
    ids[last] = max;
    bin.bic_u16(&mut ids[1..], u32::from(len).wrapping_sub(2), min, max)?;
    Ok(n)
}

#[cfg(test)]
#[path = "bitmagic_tests.rs"]
mod tests;
