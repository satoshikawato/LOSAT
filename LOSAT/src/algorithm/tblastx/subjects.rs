//! The subjects of a TBLASTX search: the records of `-subject`, or the stored sequences of a
//! `-db` database (DB-5, dbsearch-20261012 DESIGN_db5_tblastx_db.md).

use std::cell::Cell;

use crate::blastinput::fasta_reader::FastaRecord;
use crate::sequence::packed_subject::PackedNucleotide;

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
    if subject
        .iter()
        .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return None;
    }
    let ncbi4na: Vec<u8> = subject
        .iter()
        .map(|&base| match base.to_ascii_uppercase() {
            b'U' => 8,
            upper => crate::core::blast_encoding::iupacna_to_ncbi4na(upper).unwrap_or(0),
        })
        .collect();
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
