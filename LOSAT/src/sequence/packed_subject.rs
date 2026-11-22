//! A nucleotide subject held as NCBI's packed ncbi2na with its ambiguity runs and lowercase
//! intervals, in place of its letters (followups-20261014 DESIGN_subject_store.md).
//!
//! The shape is that of a BLAST database volume (`seqdb::volume::SequenceBlock`:
//! `nucleotide_packed`, `ambiguity_runs`), so a `-db` search can fill the same type. A `-subject`
//! record keeps about a quarter of a byte per base instead of one, and every form the programs
//! read is made again from it exactly: the record's letters, their ncbi4na and blastna codes, and
//! the packed bytes the seed stages read.

use std::sync::OnceLock;

use crate::core::blast_encoding::{
    encode_iupac_to_blastna, iupacna_to_ncbi4na, resolve_ncbi4na_base, NcbiRandom,
};
use crate::seqdb::volume::AmbiguityRun;
use crate::utils::dust::MaskedInterval;

/// ncbi4na code -> IUPAC letter (NCBI's `NCBI4NA_TO_IUPACNA`, upper case).
const NCBI4NA_TO_IUPACNA: [u8; 16] = *b"-ACMGRSVTWYHKDBN";

/// `letter_codes()` of a byte without an ncbi4na code.
const NO_CODE: u8 = 0xFF;
/// The bit of `letter_codes()` set for a lower-case letter.
const LOWER: u8 = 0x10;

/// Byte -> its ncbi4na code (`iupacna_to_ncbi4na`), with `LOWER` for a lower-case letter, or
/// `NO_CODE`.
fn letter_codes() -> &'static [u8; 256] {
    static CODES: OnceLock<[u8; 256]> = OnceLock::new();
    CODES.get_or_init(|| {
        std::array::from_fn(|byte| {
            let letter = byte as u8;
            match iupacna_to_ncbi4na(letter) {
                Some(code) if letter.is_ascii_lowercase() => code | LOWER,
                Some(code) => code,
                None => NO_CODE,
            }
        })
    })
}

/// Letter -> ncbi2na base for `A`, `C`, `G`, `T` (`UPPER_BASES`) or `a`, `c`, `g`, `t`
/// (`LOWER_BASES`); `AMBIGUOUS` for every other byte.
const UPPER_BASES: [u8; 256] = bases_of(*b"ACGT");
const LOWER_BASES: [u8; 256] = bases_of(*b"acgt");

const fn bases_of(letters: [u8; 4]) -> [u8; 256] {
    let mut bases = [AMBIGUOUS; 256];
    let mut base = 0;
    while base < 4 {
        bases[letters[base] as usize] = base as u8;
        base += 1;
    }
    bases
}

// No NCBI counterpart: LOSAT's speed (`bases_into`); the letters of the four 2-bit codes.
/// Packed byte -> its four bases as letters, the first base's (the high bits) first.
const BASES_OF_BYTE: [[u8; 4]; 256] = {
    let mut table = [[0u8; 4]; 256];
    let mut byte = 0;
    while byte < 256 {
        let mut slot = 0;
        while slot < 4 {
            table[byte][slot] = b"ACGT"[(byte >> (6 - 2 * slot)) & 3];
            slot += 1;
        }
        byte += 1;
    }
    table
};

/// The letters of one block of `PackedNucleotide::from_subject_letters`'s fast path.
const BLOCK: usize = 32;

/// ncbi4na code -> its blastna code (`encode_iupac_to_blastna` of the code's IUPAC letter).
fn blastna_of_ncbi4na() -> [u8; 16] {
    std::array::from_fn(|code| encode_iupac_to_blastna(&[NCBI4NA_TO_IUPACNA[code]])[0])
}

// No NCBI counterpart: LOSAT's speed (`blastna_into_par`); the packed bytes of one piece, 128 Ki
// bases (a piece is a few tens of microseconds of decoding; 3 Mb is 23 pieces).
/// The packed bytes that `blastna_into_par` decodes as one piece.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const DECODE_PIECE: usize = 32 * 1024;

// No NCBI counterpart: LOSAT's speed (the check of a block's letters is vectorized, and four
// letters take one multiplication instead of four table reads); the bytes are those of the
// letter-by-letter packing.
/// The packed bytes of the `BLOCK` letters of `block` into `out` (`BLOCK / 4` bytes) when every
/// letter is `A`, `C`, `G` or `T`, all upper case or (`lower`) all lower case; false, with `out`
/// unchanged, otherwise.
#[inline]
fn pack_block(block: &[u8], lower: bool, out: &mut [u8]) -> bool {
    let [a, c, g, t] = if lower { *b"acgt" } else { *b"ACGT" };
    let other = block.iter().fold(0u8, |other, &letter| {
        other | u8::from((letter != a) & (letter != c) & (letter != g) & (letter != t))
    });
    if other != 0 {
        return false;
    }
    // The 2-bit code of A, C, G or T in either case is bits 1-2 of the letter xor its bits 2-3
    // (0x41 0, 0x43 1, 0x47 2, 0x54 3; 0x61 0, 0x63 1, 0x67 2, 0x74 3). The multiplication moves
    // the four codes of a little-endian word to its high byte, the first letter's in the high bits
    // (the other products land below bit 24, or past bit 31 where the wrapping multiplication drops
    // them, without overlapping).
    for (byte, four) in out.iter_mut().zip(block.chunks_exact(4)) {
        let word = u32::from_le_bytes([four[0], four[1], four[2], four[3]]);
        let codes = ((word >> 1) ^ (word >> 2)) & 0x0303_0303;
        *byte = (codes.wrapping_mul(0x4010_0401) >> 24) as u8;
    }
    true
}

/// `NCBI4NA_TO_BASE` of an ambiguous code.
const AMBIGUOUS: u8 = 0xFF;
/// ncbi4na code -> ncbi2na base for `A`, `C`, `G`, `T`; `AMBIGUOUS` for the others.
const NCBI4NA_TO_BASE: [u8; 16] = {
    let mut bases = [AMBIGUOUS; 16];
    bases[1] = 0;
    bases[2] = 1;
    bases[4] = 2;
    bases[8] = 3;
    bases
};

/// One nucleotide subject: ncbi2na packed as NCBI packs a subject (`length / 4 + 1` bytes, 4
/// bases per byte from the high bits, the remainder count in the low 2 bits of the last byte),
/// the runs of letters that are not `A`, `C`, `G` or `T` (as ncbi4na codes), and the intervals
/// of lower-case letters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
/// The packed bases of a `-subject` record are those NCBI's subject set-up makes: an ambiguous
/// base is resolved by a `CRandom` made with `length`, in order (`encode_subject_ncbi2na_packed`).
/// The runs keep the real codes, which the traceback and the reports read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedNucleotide {
    packed: Box<[u8]>,
    length: usize,
    ambiguity: Box<[AmbiguityRun]>,
    lowercase: Box<[MaskedInterval]>,
}

impl PackedNucleotide {
    /// The store of a record's letters (IUPAC nucleotide letters and `-`, in either case), or
    /// `None` when a byte has no ncbi4na code (only records built by an API caller can hold one;
    /// such a search keeps the letters) or the record is longer than a run offset can hold.
    ///
    /// The letters in order, 32 or 4 at a time: the packed bytes are those of
    /// `encode_subject_ncbi2na_packed` (the random draws of the ambiguous bases in the same order),
    /// made without its two intermediate arrays of one byte per base.
    pub fn from_subject_letters(seq: &[u8]) -> Option<Self> {
        let length = seq.len();
        u32::try_from(length).ok()?;
        let letter_codes = letter_codes();
        let mut packed = vec![0u8; if length == 0 { 0 } else { length / 4 + 1 }];
        let mut ambiguity: Vec<AmbiguityRun> = Vec::new();
        let mut lowercase: Vec<MaskedInterval> = Vec::new();
        let mut lower_start: Option<usize> = None;
        // drawn from only by ambiguous bases, so made at the first one
        let mut random: Option<NcbiRandom> = None;
        let groups = length.div_ceil(4);
        let mut group = 0;
        // the group before which no block is tried (one tried block per `BLOCK` letters that are
        // not all A, C, G, T in one case, such as a run of N)
        let mut no_block_before = 0;
        while group < groups {
            let letters = &seq[group * 4..(group * 4 + 4).min(length)];
            // `BLOCK` bases A, C, G, T in the case of the letters before them: no run, no interval
            if group >= no_block_before {
                if let Some(block) = seq.get(group * 4..group * 4 + BLOCK) {
                    if pack_block(
                        block,
                        lower_start.is_some(),
                        &mut packed[group..group + BLOCK / 4],
                    ) {
                        group += BLOCK / 4;
                        continue;
                    }
                    no_block_before = group + BLOCK / 4;
                }
            }
            // four bases A, C, G, T in the case of the letters before them: no run, no interval
            if let [first, second, third, fourth] = *letters {
                let bases = if lower_start.is_none() {
                    &UPPER_BASES
                } else {
                    &LOWER_BASES
                };
                let [first, second, third, fourth] =
                    [first, second, third, fourth].map(|letter| bases[usize::from(letter)]);
                if first | second | third | fourth < 4 {
                    packed[group] = first << 6 | second << 4 | third << 2 | fourth;
                    group += 1;
                    continue;
                }
            }
            let mut byte = 0u8;
            for (slot, &letter) in letters.iter().enumerate() {
                let entry = letter_codes[usize::from(letter)];
                if entry == NO_CODE {
                    return None;
                }
                let index = group * 4 + slot;
                let code = entry & 0x0F;
                let mut base = NCBI4NA_TO_BASE[usize::from(code)];
                if base == AMBIGUOUS {
                    base = resolve_ncbi4na_base(
                        code,
                        random.get_or_insert_with(|| NcbiRandom::new(length as u32)),
                    );
                    match ambiguity.last_mut() {
                        Some(run)
                            if run.ncbi4na == code
                                && (run.offset + run.length) as usize == index =>
                        {
                            run.length += 1
                        }
                        _ => ambiguity.push(AmbiguityRun {
                            offset: index as u32,
                            length: 1,
                            ncbi4na: code,
                        }),
                    }
                }
                byte |= base << (6 - 2 * slot);
                if (entry & LOWER != 0) != lower_start.is_some() {
                    match lower_start.take() {
                        Some(start) => lowercase.push(MaskedInterval::new(start, index)),
                        None => lower_start = Some(index),
                    }
                }
            }
            packed[group] = byte;
            group += 1;
        }
        if let Some(last) = packed.last_mut() {
            *last |= (length % 4) as u8;
        }
        if let Some(start) = lower_start {
            lowercase.push(MaskedInterval::new(start, length));
        }
        Some(Self {
            packed: packed.into_boxed_slice(),
            length,
            ambiguity: ambiguity.into_boxed_slice(),
            lowercase: lowercase.into_boxed_slice(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:834-842
    /// A sequence of a BLAST database as stored (`seqdb::volume::SequenceBlock`): `packed` is
    /// its 2-bit data (`nucleotide_packed`: the bases makeblastdb stored, each ambiguous base
    /// replaced by one of its bases, in the layout of a subject's packed bytes), `length` its
    /// number of bases (`nucleotide_len`) and `runs` its ambiguity runs in stored order
    /// (`ambiguity_runs`). No lowercase intervals.
    ///
    /// The runs are applied to the letters as NCBI's reader applies them to the whole sequence
    /// (`nucleotide_blastna`): in stored order, a run starting at or past the end stops the
    /// others, a later run overwrites an earlier one. Runs in order and apart are kept as they
    /// are (makeblastdb writes them so); otherwise the runs are made again from the letters
    /// that applying them gives: the bases whose ncbi4na code differs from the stored base's.
    /// So the letters, ncbi4na and blastna codes are those of NCBI's traceback, and the packed
    /// bytes those of its preliminary stage. A run may hold the code of `A`, `C`, `G` or `T`
    /// (another letter than the stored base); it is kept, as the letters show it.
    #[cfg(test)]
    pub fn from_stored(packed: &[u8], length: usize, runs: &[AmbiguityRun]) -> Self {
        Self::from_stored_box(packed.into(), length, runs)
    }

    /// `from_stored` of the packed bytes held in `packed`, a
    /// buffer of a database read cut to them (`seqdb::SequenceBlock::x_into_nucleotide_packed`):
    /// the store keeps the buffer instead of a copy (glibc shrinks a mapped buffer in place, and
    /// the pages beyond the packed bytes go back). No NCBI counterpart: only where the bytes are
    /// held.
    pub fn x_from_stored_buffer(packed: Vec<u8>, length: usize, runs: &[AmbiguityRun]) -> Self {
        Self::from_stored_box(packed.into_boxed_slice(), length, runs)
    }

    /// `from_stored` of packed bytes already held in `packed`.
    fn from_stored_box(packed: Box<[u8]>, length: usize, runs: &[AmbiguityRun]) -> Self {
        let mut ambiguity = runs.to_vec();
        x_stored_runs(&packed, length, &mut ambiguity);
        Self {
            packed,
            length,
            ambiguity: ambiguity.into_boxed_slice(),
            lowercase: Box::default(),
        }
    }

    /// The packed bases and the ambiguity runs, borrowed (`PackedRef`): the readers of letters
    /// and bases below are the view's.
    #[inline]
    pub fn view(&self) -> PackedRef<'_> {
        PackedRef {
            packed: &self.packed,
            length: self.length,
            ambiguity: &self.ambiguity,
        }
    }

    /// The number of bases.
    pub fn len(&self) -> usize {
        self.length
    }

    /// True for a record without letters.
    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// The packed ncbi2na bytes (empty for a record without letters).
    pub fn packed(&self) -> &[u8] {
        &self.packed
    }

    /// The runs of letters other than `A`, `C`, `G`, `T`, in order, each of one ncbi4na code.
    pub fn ambiguity_runs(&self) -> &[AmbiguityRun] {
        &self.ambiguity
    }

    /// The intervals `[start, end)` of lower-case letters, in order.
    pub fn lowercase(&self) -> &[MaskedInterval] {
        &self.lowercase
    }

    /// The bytes held (packed bases, runs, intervals), for the memory accounting of a caller.
    pub fn heap_bytes(&self) -> usize {
        self.packed.len()
            + std::mem::size_of_val(&*self.ambiguity)
            + std::mem::size_of_val(&*self.lowercase)
    }

    /// Clears `out` and writes one byte per base from the packed codes through `map`.
    fn unpack_into(&self, out: &mut Vec<u8>, map: [u8; 4]) {
        self.view().unpack_into(out, map)
    }

    /// Overwrites the bases of every ambiguity run in `out` with `map[ncbi4na]`.
    fn apply_runs(&self, out: &mut [u8], map: &[u8; 16]) {
        for run in self.ambiguity.iter() {
            let start = run.offset as usize;
            out[start..start + run.length as usize].fill(map[usize::from(run.ncbi4na)]);
        }
    }

    /// The record's letters exactly as the reader gave them: IUPAC letters (and `-`) from the
    /// codes, lower case inside the lower-case intervals.
    pub fn ascii_into(&self, out: &mut Vec<u8>) {
        self.unpack_into(out, *b"ACGT");
        self.apply_runs(out, &NCBI4NA_TO_IUPACNA);
        for interval in self.lowercase.iter() {
            out[interval.start..interval.end].make_ascii_lowercase();
        }
    }

    /// The letters `[from, to)` in upper case (the letters a report row shows).
    #[inline]
    pub fn letters_into(&self, from: usize, to: usize, out: &mut Vec<u8>) {
        self.view().letters_into(from, to, out)
    }

    /// The 2-bit codes (0 `A`, 1 `C`, 2 `G`, 3 `T`; an ambiguous base as resolved) of the bases
    /// `[from, from + out.len())`.
    #[cfg(test)]
    pub fn codes_into(&self, from: usize, out: &mut [u8]) {
        assert!(from + out.len() <= self.length, "codes out of the subject");
        // the bases before the first whole byte, the whole bytes, the bases after them
        let to = from + out.len();
        let head = ((4 - from % 4) % 4).min(out.len());
        let (head_out, rest) = out.split_at_mut(head);
        for (slot, code) in head_out.iter_mut().enumerate() {
            *code = self.view().code(from + slot);
        }
        let whole = from + head;
        let mut bytes = rest.chunks_exact_mut(4);
        for (four, &byte) in (&mut bytes).zip(&self.packed[whole / 4..]) {
            four.copy_from_slice(&[byte >> 6, (byte >> 4) & 3, (byte >> 2) & 3, byte & 3]);
        }
        let tail = bytes.into_remainder();
        let tail_from = to - tail.len();
        for (slot, code) in tail.iter_mut().enumerate() {
            *code = self.view().code(tail_from + slot);
        }
    }

    /// The bases `[from, from + out.len())` as packed, one `A`, `C`, `G` or `T` each (an ambiguous
    /// base as resolved): the letters of `codes_into`'s codes.
    #[inline]
    pub fn bases_into(&self, from: usize, out: &mut [u8]) {
        self.view().bases_into(from, out)
    }

    /// The ncbi4na code of every base (`iupacna_to_ncbi4na` of each letter).
    #[cfg(test)]
    pub fn ncbi4na_into(&self, out: &mut Vec<u8>) {
        self.unpack_into(out, [1, 2, 4, 8]);
        let codes: [u8; 16] = std::array::from_fn(|code| code as u8);
        self.apply_runs(out, &codes);
    }

    /// The blastna code of every base (`encode_iupac_to_blastna` of the letters).
    pub fn blastna_into(&self, out: &mut Vec<u8>) {
        self.unpack_into(out, [0, 1, 2, 3]);
        self.apply_runs(out, &blastna_of_ncbi4na());
    }

    // No NCBI counterpart: LOSAT's speed. The traceback of a long subject decoded it on one thread
    // while the other threads of the pool waited (blastn `-subject`, 31,168 cDNAs x 3 Mb: 91 query
    // batches, each decoding the 3 Mb subject again). The pieces write disjoint bytes; the codes
    // are those of `blastna_into`.
    /// `blastna_into` with the whole packed bytes decoded in pieces of `DECODE_PIECE` bytes on the
    /// rayon pool this thread runs in (a subject of one piece, up to 128 Ki bases, is decoded on
    /// this thread). Call it inside the search pool (`SearchPool::install`).
    pub fn blastna_into_par(&self, out: &mut Vec<u8>) {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        {
            use rayon::prelude::*;
            let whole = self.length / 4;
            out.clear();
            out.reserve(self.length);
            let spare = &mut out.spare_capacity_mut()[..self.length];
            let (head, tail) = spare.split_at_mut(whole * 4);
            head.par_chunks_mut(4 * DECODE_PIECE)
                .zip(self.packed[..whole].par_chunks(DECODE_PIECE))
                .for_each(|(bases, packed)| {
                    for (four, &byte) in bases.chunks_exact_mut(4).zip(packed) {
                        for (slot, base) in four.iter_mut().enumerate() {
                            base.write((byte >> (6 - 2 * slot)) & 3);
                        }
                    }
                });
            for (slot, base) in tail.iter_mut().enumerate() {
                base.write(self.view().code(whole * 4 + slot));
            }
            // SAFETY: the pieces wrote bytes `[0, whole * 4)` and the loop above
            // `[whole * 4, self.length)`, within the capacity reserved above.
            unsafe { out.set_len(self.length) };
            self.apply_runs(out, &blastna_of_ncbi4na());
        }
        #[cfg(not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )))]
        self.blastna_into(out);
    }

    /// The bases as packed, one `A`, `C`, `G` or `T` per base: every ambiguous base resolved as
    /// NCBI's subject set-up resolves it (the draws the live TBLASTN subject translation makes
    /// for a record's letters, `x_subject_ncbi2na`).
    #[cfg(test)]
    pub fn resolved_letters_into(&self, out: &mut Vec<u8>) {
        self.unpack_into(out, *b"ACGT");
    }
}

/// A packed nucleotide sequence borrowed where it is held: a store's bytes and runs
/// (`PackedNucleotide::view`), or a sequence of a database block, read in place
/// (`PackedRef::x_stored`). It gives the letters and bases that the store gives.
#[derive(Debug, Clone, Copy)]
pub struct PackedRef<'a> {
    packed: &'a [u8],
    length: usize,
    ambiguity: &'a [AmbiguityRun],
}

impl<'a> PackedRef<'a> {
    /// A sequence of a BLAST database as `PackedNucleotide::from_stored` keeps it, without the
    /// store: `packed` its 2-bit data, `length` its number of bases and `ambiguity` the runs of
    /// `x_stored_runs`. No NCBI counterpart: only where the bytes are held.
    pub fn x_stored(packed: &'a [u8], length: usize, ambiguity: &'a [AmbiguityRun]) -> Self {
        Self {
            packed,
            length,
            ambiguity,
        }
    }

    /// The number of bases.
    pub fn len(&self) -> usize {
        self.length
    }

    /// The runs of letters other than `A`, `C`, `G`, `T`, in order, each of one ncbi4na code.
    pub fn ambiguity_runs(&self) -> &'a [AmbiguityRun] {
        self.ambiguity
    }

    /// The 2-bit code (0-3) of base `index` as packed (an ambiguous base as resolved).
    #[inline]
    fn code(&self, index: usize) -> u8 {
        (self.packed[index >> 2] >> (6 - 2 * (index & 3))) & 3
    }

    /// Clears `out` and writes one byte per base from the packed codes through `map`.
    fn unpack_into(&self, out: &mut Vec<u8>, map: [u8; 4]) {
        out.clear();
        out.reserve(self.length);
        let whole = self.length / 4;
        for &byte in &self.packed[..whole] {
            out.extend_from_slice(&[
                map[usize::from(byte >> 6)],
                map[usize::from((byte >> 4) & 3)],
                map[usize::from((byte >> 2) & 3)],
                map[usize::from(byte & 3)],
            ]);
        }
        for index in whole * 4..self.length {
            out.push(map[usize::from(self.code(index))]);
        }
    }

    /// The letters `[from, to)` in upper case (the letters a report row shows).
    pub fn letters_into(&self, from: usize, to: usize, out: &mut Vec<u8>) {
        assert!(
            from <= to && to <= self.length,
            "letters out of the subject"
        );
        out.clear();
        out.resize(to - from, 0);
        self.letters_to(from, out);
    }

    // No NCBI counterpart: LOSAT's speed (no buffer of its own); the letters of `letters_into`.
    /// `letters_into` of the letters `[from, from + out.len())`, into `out`.
    pub fn letters_to(&self, from: usize, out: &mut [u8]) {
        let to = from + out.len();
        self.bases_into(from, out);
        // the runs are in order and do not overlap: the first that ends after `from` on
        let first = self
            .ambiguity
            .partition_point(|run| (run.offset + run.length) as usize <= from);
        for run in self.ambiguity[first..].iter() {
            let start = (run.offset as usize).max(from);
            if start >= to {
                break;
            }
            let end = (run.offset as usize + run.length as usize).min(to);
            out[start - from..end - from].fill(NCBI4NA_TO_IUPACNA[usize::from(run.ncbi4na)]);
        }
    }

    // No NCBI counterpart: LOSAT's speed (one table read per packed byte); the letters are those of
    // `codes_into`'s codes.
    /// The bases `[from, from + out.len())` as packed, one `A`, `C`, `G` or `T` each (an ambiguous
    /// base as resolved): the letters of `codes_into`'s codes.
    #[inline]
    pub fn bases_into(&self, from: usize, out: &mut [u8]) {
        assert!(from + out.len() <= self.length, "bases out of the subject");
        // the bases before the first whole byte, the whole bytes, the bases after them
        let to = from + out.len();
        let head = ((4 - from % 4) % 4).min(out.len());
        let (head_out, rest) = out.split_at_mut(head);
        for (slot, base) in head_out.iter_mut().enumerate() {
            *base = b"ACGT"[usize::from(self.code(from + slot))];
        }
        let whole = from + head;
        let mut bytes = rest.chunks_exact_mut(4);
        for (four, &byte) in (&mut bytes).zip(&self.packed[whole / 4..]) {
            four.copy_from_slice(&BASES_OF_BYTE[usize::from(byte)]);
        }
        let tail = bytes.into_remainder();
        let tail_from = to - tail.len();
        for (slot, base) in tail.iter_mut().enumerate() {
            *base = b"ACGT"[usize::from(self.code(tail_from + slot))];
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:834-842 (as
// `PackedNucleotide::from_stored`): NCBI applies the runs to the whole sequence it reads.
/// The ambiguity runs of a stored sequence as `PackedNucleotide::from_stored` keeps them, made
/// in place from `runs` (its runs in stored order; `packed` its 2-bit data): clipped to the
/// `length` bases with their codes cut to 4 bits, and made again from the letters when they are
/// not in order and apart. A store and a view of a block (`PackedRef::x_stored`) read the same
/// runs.
pub fn x_stored_runs(packed: &[u8], length: usize, runs: &mut Vec<AmbiguityRun>) {
    let mut kept = 0;
    let mut in_order = true;
    let mut previous_end = 0usize;
    for index in 0..runs.len() {
        let run = runs[index];
        let position = i64::from(run.offset as i32);
        let row_end = position + i64::from(run.length);
        if row_end <= 0 {
            continue;
        }
        if position >= length as i64 {
            break;
        }
        let start = position.max(0) as usize;
        let end = row_end.min(length as i64) as usize;
        in_order &= start >= previous_end;
        previous_end = previous_end.max(end);
        runs[kept] = AmbiguityRun {
            offset: start as u32,
            length: (end - start) as u32,
            ncbi4na: run.ncbi4na & 0x0F,
        };
        kept += 1;
    }
    runs.truncate(kept);
    if in_order {
        return;
    }
    let mut codes = Vec::new();
    PackedRef::x_stored(packed, length, &[]).unpack_into(&mut codes, [1, 2, 4, 8]);
    let stored = codes.clone();
    for run in runs.iter() {
        let start = run.offset as usize;
        codes[start..start + run.length as usize].fill(run.ncbi4na);
    }
    runs.clear();
    for (index, (&code, &base)) in codes.iter().zip(&stored).enumerate() {
        if code == base {
            continue;
        }
        match runs.last_mut() {
            Some(run) if run.ncbi4na == code && (run.offset + run.length) as usize == index => {
                run.length += 1
            }
            _ => runs.push(AmbiguityRun {
                offset: index as u32,
                length: 1,
                ncbi4na: code,
            }),
        }
    }
}

/// The nucleotides of one subject as a reader of windows sees them: the record's letters, or its
/// packed store, from which a window's letters are made when it is read. A reader that ignores
/// the case of a letter (a codon's residue, a complement) reads the same bases from both.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubjectNt<'a> {
    /// The record's letters.
    Letters(&'a [u8]),
    /// The record's store.
    Packed(&'a PackedNucleotide),
}

impl<'a> From<&'a [u8]> for SubjectNt<'a> {
    fn from(letters: &'a [u8]) -> Self {
        Self::Letters(letters)
    }
}

impl<'a> From<&'a Vec<u8>> for SubjectNt<'a> {
    fn from(letters: &'a Vec<u8>) -> Self {
        Self::Letters(letters)
    }
}

impl<'a> SubjectNt<'a> {
    /// The number of bases.
    pub fn len(&self) -> usize {
        match self {
            Self::Letters(letters) => letters.len(),
            Self::Packed(store) => store.len(),
        }
    }

    /// The letters `[from, to)`: the record's own, or made in upper case in `buffer`.
    pub fn fetch<'b>(&self, from: usize, to: usize, buffer: &'b mut Vec<u8>) -> &'b [u8]
    where
        'a: 'b,
    {
        match self {
            Self::Letters(letters) => &letters[from..to],
            Self::Packed(store) => {
                store.letters_into(from, to, buffer);
                buffer
            }
        }
    }
}

#[cfg(test)]
#[path = "packed_subject_tests.rs"]
mod tests;
