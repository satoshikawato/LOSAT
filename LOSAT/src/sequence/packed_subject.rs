//! A nucleotide subject held as NCBI's packed ncbi2na with its ambiguity runs and lowercase
//! intervals, in place of its letters (followups-20261014 DESIGN_subject_store.md).
//!
//! The shape is that of a BLAST database volume (`seqdb::volume::SequenceBlock`:
//! `nucleotide_packed`, `ambiguity_runs`), so a `-db` search can fill the same type. A `-subject`
//! record keeps about a quarter of a byte per base instead of one, and every form the programs
//! read is made again from it exactly: the record's letters, their ncbi4na and blastna codes, and
//! the packed bytes the seed stages read.

use crate::core::blast_encoding::{encode_iupac_to_blastna, iupacna_to_ncbi4na};
use crate::seqdb::volume::AmbiguityRun;
use crate::utils::dust::MaskedInterval;

/// ncbi4na code -> IUPAC letter (NCBI's `NCBI4NA_TO_IUPACNA`, upper case).
const NCBI4NA_TO_IUPACNA: [u8; 16] = *b"-ACMGRSVTWYHKDBN";

/// One nucleotide subject: ncbi2na packed as NCBI packs a subject (`length / 4 + 1` bytes, 4
/// bases per byte from the high bits, the remainder count in the low 2 bits of the last byte),
/// the runs of letters that are not `A`, `C`, `G` or `T` (as ncbi4na codes), and the intervals
/// of lower-case letters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
/// ```c
///     virtual SBlastSequence GetCompressedPlusStrand() {
///         SBlastSequence retval(size());
///         string ncbi4na = kEmptyStr;
///         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
///         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
///         return retval;
/// ```
/// The packed bases of a `-subject` record are those NCBI's subject set-up makes: an ambiguous
/// base is resolved by `CRandom(length)` in order (`encode_subject_ncbi2na_packed`). The runs keep
/// the real codes, which the traceback and the reports read.
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
    pub fn from_subject_letters(seq: &[u8]) -> Option<Self> {
        u32::try_from(seq.len()).ok()?;
        let mut ambiguity: Vec<AmbiguityRun> = Vec::new();
        let mut lowercase: Vec<MaskedInterval> = Vec::new();
        let mut lower_start: Option<usize> = None;
        for (index, &letter) in seq.iter().enumerate() {
            let code = iupacna_to_ncbi4na(letter)?;
            if !matches!(code, 1 | 2 | 4 | 8) {
                match ambiguity.last_mut() {
                    Some(run)
                        if run.ncbi4na == code && (run.offset + run.length) as usize == index =>
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
            match (lower_start, letter.is_ascii_lowercase()) {
                (None, true) => lower_start = Some(index),
                (Some(start), false) => {
                    lowercase.push(MaskedInterval::new(start, index));
                    lower_start = None;
                }
                _ => {}
            }
        }
        if let Some(start) = lower_start {
            lowercase.push(MaskedInterval::new(start, seq.len()));
        }
        Some(Self {
            packed: crate::core::blast_encoding::encode_subject_ncbi2na_packed(seq)
                .into_boxed_slice(),
            length: seq.len(),
            ambiguity: ambiguity.into_boxed_slice(),
            lowercase: lowercase.into_boxed_slice(),
        })
    }

    /// The number of bases.
    pub fn len(&self) -> usize {
        self.length
    }

    /// True for a record without letters.
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
    pub fn letters_into(&self, from: usize, to: usize, out: &mut Vec<u8>) {
        assert!(
            from <= to && to <= self.length,
            "letters out of the subject"
        );
        out.clear();
        out.extend((from..to).map(|index| b"ACGT"[usize::from(self.code(index))]));
        for run in self.ambiguity.iter() {
            let start = (run.offset as usize).max(from);
            let end = (run.offset as usize + run.length as usize).min(to);
            if start < end {
                out[start - from..end - from].fill(NCBI4NA_TO_IUPACNA[usize::from(run.ncbi4na)]);
            }
        }
    }

    /// The ncbi4na code of every base (`iupacna_to_ncbi4na` of each letter).
    pub fn ncbi4na_into(&self, out: &mut Vec<u8>) {
        self.unpack_into(out, [1, 2, 4, 8]);
        let codes: [u8; 16] = std::array::from_fn(|code| code as u8);
        self.apply_runs(out, &codes);
    }

    /// The blastna code of every base (`encode_iupac_to_blastna` of the letters).
    pub fn blastna_into(&self, out: &mut Vec<u8>) {
        self.unpack_into(out, [0, 1, 2, 3]);
        let blastna: [u8; 16] =
            std::array::from_fn(|code| encode_iupac_to_blastna(&[NCBI4NA_TO_IUPACNA[code]])[0]);
        self.apply_runs(out, &blastna);
    }

    /// The bases as packed, one `A`, `C`, `G` or `T` per base: every ambiguous base resolved as
    /// NCBI's subject set-up resolves it (the letters TBLASTN's preliminary stage translates,
    /// `resolve_local_subject_ncbi2na`).
    pub fn resolved_letters_into(&self, out: &mut Vec<u8>) {
        self.unpack_into(out, *b"ACGT");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every byte the reader can give a nucleotide record: IUPAC letters and `-`, both cases.
    const LETTERS: &[u8] = b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn";

    fn random_letters(length: usize, seed: u64, ambiguous_share: u64) -> Vec<u8> {
        let mut state = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 33
        };
        let mut out = Vec::with_capacity(length);
        while out.len() < length {
            let letter = if next() % 100 < ambiguous_share {
                LETTERS[(next() % LETTERS.len() as u64) as usize]
            } else {
                b"ACGTacgt"[(next() % 8) as usize]
            };
            // runs of one letter, so that the runs and intervals are longer than one base
            let repeat = 1 + (next() % 5) as usize;
            for _ in 0..repeat.min(length - out.len()) {
                out.push(letter);
            }
        }
        out
    }

    fn cases() -> Vec<Vec<u8>> {
        let mut cases: Vec<Vec<u8>> = vec![
            Vec::new(),
            b"A".to_vec(),
            b"n".to_vec(),
            b"AC".to_vec(),
            b"ACG".to_vec(),
            b"ACGT".to_vec(),
            b"ACGTN".to_vec(),
            b"NNNNACGTNNNN".to_vec(),
            b"acgtACGTnnNN--RYry".to_vec(),
            LETTERS.to_vec(),
        ];
        // every letter at both ends and alone, with lengths across the remainder byte
        for &letter in LETTERS {
            for length in 1..=9 {
                let mut seq = vec![b'A'; length];
                seq[0] = letter;
                cases.push(seq.clone());
                seq[length - 1] = letter;
                cases.push(seq);
            }
        }
        for (seed, share) in [(1, 0), (2, 5), (3, 50), (4, 100)] {
            for length in [1, 2, 3, 4, 5, 63, 64, 65, 1000, 4099] {
                cases.push(random_letters(length, seed * 10_000 + length as u64, share));
            }
        }
        cases
    }

    #[test]
    fn every_form_is_the_one_made_from_the_letters() {
        for seq in cases() {
            let store = PackedNucleotide::from_subject_letters(&seq).expect("IUPAC letters");
            assert_eq!(store.len(), seq.len());
            assert_eq!(store.is_empty(), seq.is_empty());
            let mut out = Vec::new();
            store.ascii_into(&mut out);
            assert_eq!(out, seq, "letters of {:?}", String::from_utf8_lossy(&seq));
            assert_eq!(
                store.packed(),
                crate::core::blast_encoding::encode_subject_ncbi2na_packed(&seq).as_slice()
            );
            store.ncbi4na_into(&mut out);
            let ncbi4na: Vec<u8> = seq
                .iter()
                .map(|&letter| iupacna_to_ncbi4na(letter).unwrap())
                .collect();
            assert_eq!(out, ncbi4na);
            store.blastna_into(&mut out);
            assert_eq!(out, encode_iupac_to_blastna(&seq));
            // TBLASTN's seed input (`resolve_local_subject_ncbi2na`): the letters' ncbi4na codes
            // resolved by `CRandom(length)`, as letters
            store.resolved_letters_into(&mut out);
            let resolved: Vec<u8> =
                crate::core::blast_encoding::resolve_ncbi4na_to_ncbi2na(&ncbi4na)
                    .iter()
                    .map(|&code| b"ACGT"[usize::from(code)])
                    .collect();
            assert_eq!(out, resolved);
            let upper = seq.to_ascii_uppercase();
            for (from, to) in [
                (0, seq.len()),
                (0, seq.len() / 2),
                (seq.len() / 3, seq.len()),
            ] {
                store.letters_into(from, to, &mut out);
                assert_eq!(out, &upper[from..to]);
            }
            let lower: Vec<_> =
                crate::algorithm::blastn::coordination::collect_lowercase_masks(&seq)
                    .iter()
                    .map(|interval| (interval.start, interval.end))
                    .collect();
            let held: Vec<_> = store
                .lowercase()
                .iter()
                .map(|interval| (interval.start, interval.end))
                .collect();
            assert_eq!(held, lower);
        }
    }

    #[test]
    fn runs_are_maximal_runs_of_one_code() {
        let store = PackedNucleotide::from_subject_letters(b"ANNnNRRACGT--N").unwrap();
        let runs: Vec<_> = store
            .ambiguity_runs()
            .iter()
            .map(|run| (run.offset, run.length, run.ncbi4na))
            .collect();
        assert_eq!(runs, [(1, 4, 15), (5, 2, 5), (11, 2, 0), (13, 1, 15)]);
        assert_eq!(store.heap_bytes(), 14 / 4 + 1 + 4 * 12 + 16);
    }

    #[test]
    fn a_byte_without_a_code_keeps_the_letters() {
        for seq in [&b"ACGU"[..], b"ACG*", b"AC GT", b"ACGX", &[b'A', 0xFF]] {
            assert!(
                PackedNucleotide::from_subject_letters(seq).is_none(),
                "{seq:?}"
            );
        }
    }
}
