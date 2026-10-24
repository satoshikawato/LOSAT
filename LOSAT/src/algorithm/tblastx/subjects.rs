//! The subjects of a TBLASTX search: the records of `-subject`, or the stored sequences of a
//! `-db` database (DB-5, dbsearch-20261012 DESIGN_db5_tblastx_db.md).

use std::cell::Cell;
use std::sync::OnceLock;

use crate::algorithm::tblastx::translation::{x_generate_frames_in_windows, QueryFrame};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::sequence::packed_subject::PackedNucleotide;
use crate::utils::genetic_code::GeneticCode;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1403-1405
// ```c
//     /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
//       sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
//     seq_arg.encoding = eBlastEncodingProtein;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2650-2653
// ```c
//       	seq_arg.oid = subject_blk->oid;
//       	seq_arg.encoding =
//        	  (kTranslateSubject ? eBlastEncodingNcbi4na : eBlastEncodingNucleotide);
// ```
// The preliminary stage translates a subject's 2-bit bases, the ungapped re-evaluation its
// ncbi4na codes. A record of `-subject` gives its letters and the bases its ambiguous letters
// resolve to (`preliminary_subject_bases`); a database sequence gives the letters of its stored
// codes and the bases makeblastdb stored (D2).
/// The subjects of one search, by OID (`-db`: the ordinal).
#[derive(Clone, Copy)]
pub(crate) enum Subjects<'a> {
    /// The records of `-subject` (after `-subject_loc`).
    Fasta(&'a [FastaRecord]),
    /// The sequences of a `-db` database, read whole.
    Stored(&'a [PackedNucleotide]),
    /// A `-db` database whose sequences could not be read: no subject is searched, and the
    /// first batch with a valid context raises the read error (`SubjectsUnreadable`).
    Unreadable,
}

impl<'a> Subjects<'a> {
    /// The number of subjects.
    pub(crate) fn count(&self) -> usize {
        match self {
            Self::Fasta(records) => records.len(),
            Self::Stored(stores) => stores.len(),
            Self::Unreadable => 0,
        }
    }

    /// The number of bases of subject `oid`.
    pub(crate) fn length(&self, oid: usize) -> usize {
        match self {
            Self::Fasta(records) => records[oid].seq().len(),
            Self::Stored(stores) => stores[oid].len(),
            Self::Unreadable => 0,
        }
    }

    /// The numbers of bases of the subjects, in OID order.
    pub(crate) fn lengths(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.count()).map(|oid| self.length(oid))
    }

    /// `f` of the letters of subject `oid` and, when it has ambiguous bases, of the bases its
    /// preliminary stage translates (`A`, `C`, `G` or `T` for every base). A stored sequence's
    /// letters and bases are made in buffers that the thread keeps between subjects (one above
    /// `KEPT_BUFFER` is freed after use); a nested call (a pool thread that takes another
    /// subject while it waits) gets buffers of its own.
    pub(crate) fn with_search_bases<R>(
        &self,
        oid: usize,
        f: impl FnOnce(&[u8], Option<&[u8]>) -> R,
    ) -> R {
        match self {
            Self::Fasta(records) => {
                let letters = records[oid].seq();
                let bases = preliminary_subject_bases(letters);
                f(letters, bases.as_deref())
            }
            Self::Stored(stores) => {
                let store = &stores[oid];
                with_buffer(
                    &LETTERS,
                    |letters| store.ascii_into(letters),
                    |letters| {
                        if store.ambiguity_runs().is_empty() {
                            f(letters, None)
                        } else {
                            with_buffer(
                                &BASES,
                                |bases| store.resolved_letters_into(bases),
                                |bases| f(letters, Some(bases)),
                            )
                        }
                    },
                )
            }
            Self::Unreadable => unreachable!("an unreadable database has no subject to search"),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
    // ```c
    //             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
    //                                      backup.full_range.right,
    //                                      subject->gen_code_string, &translation_buffer,
    //                                      &frame_offsets, NULL);
    // ```
    // EXPERIMENT (LOSAT_X_TBXFRAMES): the frames of the bases `with_search_bases` gives the
    // preliminary stage (its `bases`, or its letters when the subject has no ambiguous base),
    // made in windows (`x_generate_frames_in_windows`) without the subject's letters, its bases
    // or their reverse complement in a buffer of the subject's length: a record's ambiguous
    // letters are resolved window by window in letter order from one `CRandom` seeded with the
    // record's length, which draws what `preliminary_subject_bases` draws; a stored sequence
    // gives the bases of its 2-bit codes (`resolved_letters_into`).
    /// The frames the preliminary stage of subject `oid` translates, and whether they are the
    /// frames of its letters as well (true when it has no ambiguous base).
    pub(crate) fn x_preliminary_frames(
        &self,
        oid: usize,
        table: &GeneticCode,
    ) -> (Vec<QueryFrame>, bool) {
        match self {
            Self::Fasta(records) => {
                let letters = records[oid].seq();
                let unambiguous = is_unambiguous(letters);
                let by_table = x_resolve_by_table();
                // A record of A, C, G and T only copies its letters with either side.
                if by_table {
                    crate::utils::xreach::x_reach!("LOSAT_X_TBXRESOLVE", !unambiguous);
                }
                if unambiguous {
                    let frames = x_generate_frames_in_windows(
                        letters.len(),
                        &mut |from, out| out.copy_from_slice(&letters[from..from + out.len()]),
                        table,
                    );
                    return (frames, true);
                }
                (x_resolved_record_frames(letters, table, by_table), false)
            }
            Self::Stored(stores) => {
                let store = &stores[oid];
                let frames = x_generate_frames_in_windows(
                    store.len(),
                    &mut |from, out| {
                        store.codes_into(from, out);
                        for base in out.iter_mut() {
                            *base = b"ACGT"[usize::from(*base)];
                        }
                    },
                    table,
                );
                (frames, store.ambiguity_runs().is_empty())
            }
            Self::Unreadable => unreachable!("an unreadable database has no subject to search"),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2650-2653
    // ```c
    //       	seq_arg.oid = subject_blk->oid;
    //       	seq_arg.encoding =
    //        	  (kTranslateSubject ? eBlastEncodingNcbi4na : eBlastEncodingNucleotide);
    // ```
    // EXPERIMENT (LOSAT_X_TBXFRAMES): the frames of the letters `with_search_bases` gives
    // (`s_frames`), made in windows: a record's letters, or a stored sequence's letters of its
    // codes (`letters_into`, the upper-case letters of `ascii_into`: a stored sequence has no
    // lower-case interval, and the translation reads both cases alike).
    /// The frames of the letters of subject `oid` (the ungapped re-evaluation's).
    pub(crate) fn x_letter_frames(&self, oid: usize, table: &GeneticCode) -> Vec<QueryFrame> {
        match self {
            Self::Fasta(records) => {
                let letters = records[oid].seq();
                x_generate_frames_in_windows(
                    letters.len(),
                    &mut |from, out| out.copy_from_slice(&letters[from..from + out.len()]),
                    table,
                )
            }
            Self::Stored(stores) => {
                let store = &stores[oid];
                let mut window: Vec<u8> = Vec::new();
                x_generate_frames_in_windows(
                    store.len(),
                    &mut |from, out| {
                        store.letters_into(from, from + out.len(), &mut window);
                        out.copy_from_slice(&window);
                    },
                    table,
                )
            }
            Self::Unreadable => unreachable!("an unreadable database has no subject to search"),
        }
    }
}

// No NCBI counterpart: the variable is read once per process; it chooses when LOSAT translates
// a subject's frames, not a value NCBI computes.
/// LOSAT_X_TBXFRAMES: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
pub(crate) fn x_frames_late() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXFRAMES"))
}

// No NCBI counterpart: the variable is read once per process; it chooses how LOSAT resolves a
// record's letters, not a value NCBI computes.
/// LOSAT_X_TBXRESOLVE: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
pub(crate) fn x_resolve_by_table() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXRESOLVE"))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:424-425,441-450
// ```c++
// static unsigned char ctable[16] = {0xFF, 0x00, 0x01, 0xFF, 0x02, 0xFF, 0xFF, 0xFF,
// 		                           0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
// ...
//         unsigned char c = ctable[b];
//
//         if (c  != 0xFF) {
//             // No ambiguities, so we can do this the easy way.
//         	ncbi2na[i] = c;
//
//         } else {
//             if (b == 0 || b == 0x0F) {
//             	//gap or N
//                 ncbi2na[i] = random.GetRand() & 0x3;
// ```
// NCBI draws from `CRandom` only for an ncbi4na value that `ctable` does not map (not 1, 2, 4 or
// 8). EXPERIMENT (LOSAT_X_TBXRESOLVE): the bytes whose value `ncbi4na_of` gives is one of those
// four read their base from `ONE_BASE`; only the others go through `ncbi4na_of` and
// `resolve_ncbi4na_base`, so the draws and the bases are those of the letter-by-letter walk.
/// The base of each byte whose ncbi4na value (`ncbi4na_of`) has one base: `A`, `C`, `G`, `T` and
/// `U` (as `T`), in either case; 0 for every other byte.
const ONE_BASE: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'A';
    table[b'a' as usize] = b'A';
    table[b'C' as usize] = b'C';
    table[b'c' as usize] = b'C';
    table[b'G' as usize] = b'G';
    table[b'g' as usize] = b'G';
    table[b'T' as usize] = b'T';
    table[b't' as usize] = b'T';
    table[b'U' as usize] = b'T';
    table[b'u' as usize] = b'T';
    table
};

/// The frames of the bases the preliminary stage translates for a record with a letter other than
/// A, C, G and T (`Subjects::x_preliminary_frames`), the letters resolved window by window from one
/// `CRandom` of the record's length; `by_table` is LOSAT_X_TBXRESOLVE.
fn x_resolved_record_frames(
    letters: &[u8],
    table: &GeneticCode,
    by_table: bool,
) -> Vec<QueryFrame> {
    let mut random = crate::core::blast_encoding::NcbiRandom::new(letters.len() as u32);
    let mut resolve = |letter: u8| {
        let code = ncbi4na_of(letter);
        b"ACGT"[usize::from(crate::core::blast_encoding::resolve_ncbi4na_base(
            code,
            &mut random,
        ))]
    };
    x_generate_frames_in_windows(
        letters.len(),
        &mut |from, out| {
            if by_table {
                for (base, &letter) in out.iter_mut().zip(&letters[from..]) {
                    let known = ONE_BASE[usize::from(letter)];
                    *base = if known != 0 { known } else { resolve(letter) };
                }
            } else {
                for (base, &letter) in out.iter_mut().zip(&letters[from..]) {
                    *base = resolve(letter);
                }
            }
        },
        table,
    )
}

/// Whether every letter is `A`, `C`, `G` or `T` in either case (`preliminary_subject_bases`
/// returns `None`).
fn is_unambiguous(letters: &[u8]) -> bool {
    letters
        .iter()
        .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

/// The ncbi4na code `preliminary_subject_bases` gives a letter.
fn ncbi4na_of(letter: u8) -> u8 {
    match letter.to_ascii_uppercase() {
        b'U' => 8,
        upper => crate::core::blast_encoding::iupacna_to_ncbi4na(upper).unwrap_or(0),
    }
}

/// The bases of a subject as the preliminary search reads them: NCBI's compressed (ncbi2na)
/// subject, in which each ambiguous base is a compatible base drawn from `CRandom` seeded
/// with the subject's length (`resolve_ncbi4na_to_ncbi2na`). `None` for a subject of
/// A, C, G and T only, whose compressed bases are its own.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
/// ```c
///     virtual SBlastSequence GetCompressedPlusStrand() {
///         SBlastSequence retval(size());
///         string ncbi4na = kEmptyStr;
///         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
///         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
///         return retval;
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:94-103
/// ```c
/// const Uint1 IUPACNA_TO_NCBI4NA[128]={
///  0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
///  0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0,
/// ```
/// The ncbi4na code of a letter is that of `GeneticCode::get` (both cases, `U` as `T`);
/// another byte is a gap (0).
pub(crate) fn preliminary_subject_bases(subject: &[u8]) -> Option<Vec<u8>> {
    if is_unambiguous(subject) {
        return None;
    }
    let ncbi4na: Vec<u8> = subject.iter().map(|&base| ncbi4na_of(base)).collect();
    Some(
        crate::core::blast_encoding::resolve_ncbi4na_to_ncbi2na(&ncbi4na)
            .into_iter()
            .map(|code| b"ACGT"[code as usize])
            .collect(),
    )
}

/// `f` of the bytes `make` writes in this thread's buffer `cell` (see `with_search_bases`).
fn with_buffer<R>(
    cell: &'static std::thread::LocalKey<Cell<Vec<u8>>>,
    make: impl FnOnce(&mut Vec<u8>),
    f: impl FnOnce(&[u8]) -> R,
) -> R {
    let mut buffer = cell.with(Cell::take);
    make(&mut buffer);
    let result = f(&buffer);
    if buffer.capacity() <= KEPT_BUFFER {
        cell.with(|kept| kept.set(buffer));
    }
    result
}

/// The largest buffer a thread keeps between subjects.
const KEPT_BUFFER: usize = 8 << 20;

thread_local! {
    /// This thread's letters of a stored subject.
    static LETTERS: Cell<Vec<u8>> = const { Cell::new(Vec::new()) };
    /// This thread's bases of a stored subject with ambiguous bases.
    static BASES: Cell<Vec<u8>> = const { Cell::new(Vec::new()) };
}

/// The error of the first batch with a valid context when the database's sequences could not be
/// read (`Subjects::Unreadable`); the caller raises the read error it kept.
#[derive(Debug)]
pub(crate) struct SubjectsUnreadable;

impl std::fmt::Display for SubjectsUnreadable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "the database's sequences could not be read")
    }
}

impl std::error::Error for SubjectsUnreadable {}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-735
// ```c
//    if (eff_len_options->db_length > 0)
//       db_length = eff_len_options->db_length;
//    else
//       db_length = eff_len_params->real_db_length;
// ...
//    if (Blast_SubjectIsTranslated(program_number))
//       db_length = db_length/3;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:909-911,974-976
// ```c
//       total_length = BlastSeqSrcGetTotLenStats(seq_src);
//       if (total_length <= 0)
//           total_length = BlastSeqSrcGetTotLen(seq_src);
// ...
//    } else {
//        min_subject_length = (Int4) (total_length/num_seqs);
//    }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407,1451-1455
// ```c
//     db_length = BlastSeqSrcGetTotLen(seq_src);
// ...
//           CalculateLinkHSPCutoffs(program_number, query_info, sbp,
//             hit_params->link_hsp_params, word_params, db_length,
//             seq_arg.seq->length);
// ```
// A subject set (seqsrc_multiseq.cpp:293-298) and a database (seqsrc_seqdb.cpp:183-187) give
// their real totals; `-dbsize` replaces only the effective lengths' database length.
/// The numbers of the subject set that the statistics read, once per run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SubjectSetNumbers {
    /// The bases of every subject (`-db`: the index's total).
    pub(crate) total: i64,
    /// The number of subjects (`-db`: the index's count).
    pub(crate) count: i64,
    /// `-dbsize`.
    pub(crate) dbsize: Option<i64>,
}

impl SubjectSetNumbers {
    /// The numbers of the records of `-subject`.
    pub(crate) fn of_records(records: &[FastaRecord], dbsize: Option<i64>) -> Self {
        Self {
            total: records.iter().map(|record| record.seq().len() as i64).sum(),
            count: records.len() as i64,
            dbsize,
        }
    }

    /// The database length of the effective lengths, in bases: `-dbsize` when above 0, else
    /// the real total.
    pub(crate) fn effective_length(&self) -> i64 {
        self.dbsize
            .filter(|&dbsize| dbsize > 0)
            .unwrap_or(self.total)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:836-843
// ```c
//         	 Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
//
//         	 // Just in case effective_db_length < 0
//         	 if (effective_db_length <= 0)
//         		 effective_db_length = 1;
//
//              effective_search_space = effective_db_length *
//                              (query_length - length_adjustment);
// ```
/// LOSAT's refusal of a search whose effective search space overflows NCBI's `Int8`
/// (`db_length` bases of the database, or of `-dbsize`, against a context of `query_length`
/// residues): NCBI's arithmetic overflows there (as BLASTP's and BLASTX's `-dbsize`).
#[derive(Debug)]
pub(crate) struct SearchSpaceOverflow {
    pub(crate) db_length: i64,
    pub(crate) query_length: usize,
}

impl std::fmt::Display for SearchSpaceOverflow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "a database length of {} letters with a query of {} letters overflows NCBI BLAST+'s 64-bit effective search space, which is not supported by LOSAT's TBLASTX",
            self.db_length, self.query_length
        )
    }
}

impl std::error::Error for SearchSpaceOverflow {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastx::translation::generate_frames;
    use crate::seqdb::volume::AmbiguityRun;

    fn assert_same(got: &[QueryFrame], expected: &[QueryFrame], what: &str) {
        assert_eq!(got.len(), expected.len(), "{what}");
        for (g, e) in got.iter().zip(expected) {
            assert_eq!(
                (g.frame, &g.aa_seq, g.aa_len, g.orig_len),
                (e.frame, &e.aa_seq, e.aa_len, e.orig_len),
                "{what}"
            );
        }
    }

    /// LOSAT_X_TBXFRAMES: the two methods against the frames of `with_search_bases` (the frames
    /// `process_subject` makes without the switch), for the subjects of both arms.
    fn check(subjects: Subjects<'_>, table: &GeneticCode) {
        for oid in 0..subjects.count() {
            let (letters, preliminary) = subjects.with_search_bases(oid, |letters, bases| {
                (
                    generate_frames(letters, table),
                    bases.map(|bases| generate_frames(bases, table)),
                )
            });
            let (frames, of_letters) = subjects.x_preliminary_frames(oid, table);
            assert_eq!(of_letters, preliminary.is_none(), "oid {oid}");
            assert_same(
                &frames,
                preliminary.as_deref().unwrap_or(&letters),
                "preliminary",
            );
            assert_same(&subjects.x_letter_frames(oid, table), &letters, "letters");
        }
    }

    fn letters(seed: u64, len: usize, alphabet: &[u8]) -> Vec<u8> {
        let mut state = seed | 1;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                alphabet[(state % alphabet.len() as u64) as usize]
            })
            .collect()
    }

    #[test]
    fn frames_of_records_are_those_of_their_letters_and_bases() {
        let mut records = Vec::new();
        for (i, len) in [0usize, 2, 3, 4, 50, 1000, 200_003].into_iter().enumerate() {
            records.push(FastaRecord::new(
                format!("a{i}"),
                b"",
                &letters(i as u64, len, b"ACGTacgt"),
            ));
            records.push(FastaRecord::new(
                format!("n{i}"),
                b"",
                &letters(i as u64 + 99, len, b"ACGTACGTacgtNNNRYKMSWBDHVnU-"),
            ));
        }
        // One N block across a window border.
        let mut block = letters(7, 400_000, b"ACGT");
        block[196_000..197_000].fill(b'N');
        records.push(FastaRecord::new("block", b"", &block));
        for id in [1, 4, 11] {
            check(Subjects::Fasta(&records), &GeneticCode::from_id(id));
        }
    }

    /// LOSAT_X_TBXRESOLVE: every cell of `ONE_BASE` against `ncbi4na_of` and `resolve_ncbi4na_base`:
    /// a byte of the table gets the resolve's base without a draw, and every other byte has an
    /// ncbi4na value of no base or of more than one.
    #[test]
    fn one_base_is_the_resolve_without_a_draw() {
        use crate::core::blast_encoding::{resolve_ncbi4na_base, NcbiRandom};
        for byte in 0..=u8::MAX {
            let code = ncbi4na_of(byte);
            let known = ONE_BASE[usize::from(byte)];
            if known == 0 {
                assert_ne!(code.count_ones(), 1, "byte {byte}");
                continue;
            }
            let (mut drawn, mut fresh) = (NcbiRandom::new(7), NcbiRandom::new(7));
            let base = b"ACGT"[usize::from(resolve_ncbi4na_base(code, &mut drawn))];
            assert_eq!(known, base, "byte {byte}");
            assert_eq!(drawn.get_rand(), fresh.get_rand(), "byte {byte}: a draw");
        }
    }

    /// LOSAT_X_TBXRESOLVE: the table walk against the letter-by-letter walk (`=0`), every frame,
    /// for records of every IUPAC letter in both cases, `U`, gaps and bytes without a code, with
    /// one ambiguous letter, and with an N block across a window border.
    #[test]
    fn the_table_walk_gives_the_letter_walks_frames() {
        let mut records = vec![
            b"N".to_vec(),
            b"acgtN".to_vec(),
            b"uUuUaN".to_vec(),
            letters(5, 3, b"ACGTN"),
        ];
        for (i, len) in [1usize, 2, 3, 4, 7, 50, 1000, 200_003]
            .into_iter()
            .enumerate()
        {
            records.push(letters(
                i as u64 + 11,
                len,
                b"ACGTacgtNnRrYyKkMmSsWwBbDdHhVvUu-*0 ",
            ));
            let mut one = letters(i as u64 + 31, len, b"ACGTacgt");
            one[len / 2] = b'R';
            records.push(one);
        }
        let mut block = letters(17, 400_000, b"ACGTacgu");
        block[196_000..197_000].fill(b'N');
        records.push(block);
        for id in [1, 4, 11] {
            let table = GeneticCode::from_id(id);
            for (i, record) in records.iter().enumerate() {
                assert_same(
                    &x_resolved_record_frames(record, &table, true),
                    &x_resolved_record_frames(record, &table, false),
                    &format!("record {i} code {id}"),
                );
            }
        }
    }

    #[test]
    fn frames_of_stored_sequences_are_those_of_their_codes_and_stored_bases() {
        let mut stores = Vec::new();
        for (i, len) in [3usize, 7, 64, 1001, 300_001].into_iter().enumerate() {
            let packed: Vec<u8> =
                letters(i as u64 + 5, len / 4 + 1, &(0..=255).collect::<Vec<u8>>());
            stores.push(PackedNucleotide::from_stored(&packed, len, &[]));
            // Runs of ambiguity codes, also of A/C/G/T codes other than the stored base, and
            // one across a window border; the stored bases under them are arbitrary (D2).
            let mut runs = vec![
                AmbiguityRun {
                    offset: 0,
                    length: 1,
                    ncbi4na: 15,
                },
                AmbiguityRun {
                    offset: 2,
                    length: 1,
                    ncbi4na: 5,
                },
            ];
            if len > 1000 {
                runs.push(AmbiguityRun {
                    offset: 500,
                    length: 3,
                    ncbi4na: 1,
                });
                runs.push(AmbiguityRun {
                    offset: 900,
                    length: 40,
                    ncbi4na: 14,
                });
            }
            if len > 200_000 {
                runs.push(AmbiguityRun {
                    offset: 196_600,
                    length: 50,
                    ncbi4na: 15,
                });
            }
            stores.push(PackedNucleotide::from_stored(&packed, len, &runs));
        }
        for id in [1, 2, 11] {
            check(Subjects::Stored(&stores), &GeneticCode::from_id(id));
        }
    }
}
