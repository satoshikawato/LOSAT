//! The subjects of a TBLASTX search: the records of `-subject`, or the stored sequences of a
//! `-db` database (DB-5, dbsearch-20261012 DESIGN_db5_tblastx_db.md).

#[cfg(test)]
use std::cell::Cell;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use anyhow::Result;

use crate::algorithm::common::x_db_stream::{XDbStream, XSequence};
use crate::algorithm::tblastx::translation::{x_generate_frames_in_windows, QueryFrame};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::seqdb::SequenceBlock;
use crate::sequence::packed_subject::{PackedNucleotide, PackedRef};
use crate::utils::genetic_code::GeneticCode;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1403-1405
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2650-2653
// The preliminary stage translates a subject's 2-bit bases, the ungapped re-evaluation its
// ncbi4na codes. A record of `-subject` gives its letters and the bases its ambiguous letters
// resolve to (`preliminary_subject_bases`); a database sequence gives the letters of its stored
// codes and the bases makeblastdb stored (D2); a record packed as read (LOSAT_X_TBXSTORE) gives
// the letters of its codes and the bases its ambiguous letters were resolved to when it was packed
// (`CRandom` of its length, as `preliminary_subject_bases`).
/// The subjects of one search, by OID (`-db`: the ordinal).
#[derive(Clone, Copy)]
pub(crate) enum Subjects<'a> {
    /// The records of `-subject` (after `-subject_loc`).
    Fasta(&'a [FastaRecord]),
    /// (LOSAT_X_TBXSTORE) The records of `-subject` (after `-subject_loc`), packed as read.
    Stored(&'a [PackedNucleotide]),
    /// The sequences of a `-db` database, read where they are searched (`XDbStream`): a short
    /// sequence's frames are made from the bytes of its block, a long one (a block of its own) is
    /// read alone and shared, and the stores are those of a copy of the database
    /// (`common::database::read_nucleotide_stores`), so every frame is the same.
    Stream(&'a XDbStream<'a>),
    /// A `-db` database whose sequences could not be read: no subject is searched, and the
    /// first batch with a valid context raises the read error (`SubjectsUnreadable`).
    Unreadable,
}

impl<'a> Subjects<'a> {
    /// The database of a `-db` search whose sequences are read (`None` for `-subject`).
    pub(crate) fn seqdb(&self) -> Option<&'a crate::seqdb::SeqDb> {
        match self {
            Self::Stream(stream) => Some(stream.seqdb()),
            Self::Fasta(_) | Self::Stored(_) | Self::Unreadable => None,
        }
    }

    /// The number of subjects.
    pub(crate) fn count(&self) -> usize {
        match self {
            Self::Fasta(records) => records.len(),
            Self::Stored(stores) => stores.len(),
            Self::Stream(stream) => stream.count(),
            Self::Unreadable => 0,
        }
    }

    /// The number of bases of subject `oid` of the subjects held in memory (a `-db` stream reads
    /// a sequence's length with the sequence, `x_preliminary_frames`).
    fn held_length(&self, oid: usize) -> usize {
        match self {
            Self::Fasta(records) => records[oid].seq().len(),
            Self::Stored(stores) => stores[oid].len(),
            Self::Stream(_) => unreachable!("a -db stream reads its lengths with its sequences"),
            Self::Unreadable => 0,
        }
    }

    /// The number of bases of the first subject (a stream's from its check); there is one.
    pub(crate) fn first_length(&self) -> usize {
        match self {
            Self::Stream(stream) => stream.first_length(),
            _ => self.held_length(0),
        }
    }

    /// The number of bases of the longest subject and of all the subjects (a stream holds them
    /// from its check).
    pub(crate) fn longest_and_total(&self) -> (usize, u64) {
        match self {
            Self::Stream(stream) => stream.longest_and_total(),
            _ => (0..self.count()).fold((0, 0), |(longest, total), oid| {
                let length = self.held_length(oid);
                (longest.max(length), total + length as u64)
            }),
        }
    }

    /// The reference of the unit tests: `f` of the letters of subject `oid` and, when it has
    /// ambiguous bases, of the bases its preliminary stage translates (`A`, `C`, `G` or `T` for
    /// every base), whole, as the search made them before its frames were made in windows
    /// (`x_preliminary_frames`, `x_letter_frames`). A stored sequence's letters and bases are made
    /// in buffers that the thread keeps between subjects (one above `KEPT_BUFFER` is freed after
    /// use); a nested call gets buffers of its own.
    #[cfg(test)]
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
            Self::Stream(stream) => {
                let store = stream.read(oid).unwrap();
                Subjects::Stored(std::slice::from_ref(&store)).with_search_bases(0, f)
            }
            Self::Unreadable => unreachable!("an unreadable database has no subject to search"),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
    // The frames of the bases `with_search_bases` gives the preliminary stage (its `bases`, or
    // its letters when the subject has no ambiguous base), made in windows
    // (`x_generate_frames_in_windows`) without the subject's letters, its bases or their reverse
    // complement in a buffer of the subject's length: a record's ambiguous letters are resolved
    // window by window in letter order from one `CRandom` seeded with the record's length, which
    // draws what `preliminary_subject_bases` draws; a stored sequence gives the bases of its
    // 2-bit codes (`bases_into`).
    /// The number of bases of subject `oid`, the frames its preliminary stage translates, whether
    /// they are the frames of its letters as well (true when it has no ambiguous base), and what
    /// the search keeps for the frames of its letters (`x_letter_frames`); none (no frame made)
    /// for a subject of fewer than `shortest` bases, which the search does not search. A `-db`
    /// stream's length is read with the sequence, in the one look-up of its block. A read of a
    /// `-db` stream that fails is an `XReadError`.
    pub(crate) fn x_preliminary_frames(
        &self,
        oid: usize,
        table: &GeneticCode,
        shortest: usize,
    ) -> Result<Option<XPreliminary>> {
        match self {
            Self::Fasta(records) if records[oid].seq().len() < shortest => Ok(None),
            Self::Stored(stores) if stores[oid].len() < shortest => Ok(None),
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
                    return Ok(Some((letters.len(), frames, true, None)));
                }
                Ok(Some((
                    letters.len(),
                    x_resolved_record_frames(letters, table, by_table),
                    false,
                    None,
                )))
            }
            Self::Stored(stores) => {
                let store = &stores[oid];
                let frames = x_generate_frames_in_windows(
                    store.len(),
                    &mut |from, out| store.bases_into(from, out),
                    table,
                );
                Ok(Some((
                    store.len(),
                    frames,
                    store.ambiguity_runs().is_empty(),
                    None,
                )))
            }
            // A short sequence keeps its block until its letters' frames are made, when it has
            // ambiguous bases (a nested subject on this thread may take another block); a long one
            // is read again then, instead of being held through the search.
            // With `-db_hard_mask` (series NP, M4b PR 2) the frames of the letters are those of
            // the stages after the preliminary one, `B` over the masks (`x_letter_frames`).
            Self::Stream(stream) => stream.with_block(oid, |sequence| {
                let length = sequence.len();
                if length < shortest {
                    return Ok(None);
                }
                let hard = stream.mask().is_some_and(|mask| mask.hard());
                let block = match &sequence {
                    XSequence::Block { block, i, .. } if !hard => Some((Arc::clone(block), *i)),
                    _ => None,
                };
                let (frames, of_letters) = sequence.with_view(|view| {
                    let frames = x_generate_frames_in_windows(
                        view.len(),
                        &mut |from, out| view.bases_into(from, out),
                        table,
                    );
                    (frames, !hard && view.ambiguity_runs().is_empty())
                });
                Ok(Some((
                    length,
                    frames,
                    of_letters,
                    block.filter(|_| !of_letters),
                )))
            }),
            Self::Unreadable => unreachable!("an unreadable database has no subject to search"),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2650-2653
    // The frames of the letters `with_search_bases` gives (`s_frames`), made in windows: a
    // record's letters, or a stored sequence's letters of its codes (`letters_to`, the
    // upper-case letters of `ascii_into`; the translation reads both cases alike).
    /// The frames of the letters of subject `oid` (the ungapped re-evaluation's); `kept` is what
    /// `x_preliminary_frames` kept of it. A read of a `-db` stream that fails is an `XReadError`.
    pub(crate) fn x_letter_frames(
        &self,
        oid: usize,
        table: &GeneticCode,
        kept: XKept,
    ) -> Result<Vec<QueryFrame>> {
        let of_view = |view: PackedRef<'_>| {
            x_generate_frames_in_windows(
                view.len(),
                &mut |from, out| view.letters_to(from, out),
                table,
            )
        };
        Ok(match self {
            Self::Fasta(records) => {
                let letters = records[oid].seq();
                x_generate_frames_in_windows(
                    letters.len(),
                    &mut |from, out| out.copy_from_slice(&letters[from..from + out.len()]),
                    table,
                )
            }
            Self::Stored(stores) => of_view(stores[oid].view()),
            Self::Stream(stream) => match kept {
                Some((block, i)) => {
                    let (length, runs) = stream.parts_in(&block, i, oid)?;
                    XSequence::Block {
                        block: &block,
                        i,
                        length,
                        runs,
                    }
                    .with_view(of_view)
                }
                // With `-db_hard_mask`, the sequence as the later stages read it.
                None => of_view(stream.read_stage(oid)?.view()),
            },
            Self::Unreadable => unreachable!("an unreadable database has no subject to search"),
        })
    }
}

/// What `Subjects::x_preliminary_frames` makes of a subject long enough to be searched: its
/// number of bases, its preliminary frames, whether they are the frames of its letters, and what
/// the search keeps for those (`XKept`).
pub(crate) type XPreliminary = (usize, Vec<QueryFrame>, bool, XKept);

/// What the search of a `-db` subject keeps from its preliminary frames to the frames of its
/// letters (`Subjects::x_letter_frames`): the block that holds a short sequence with ambiguous
/// bases, and its place in the block. `None` for every other subject.
pub(crate) type XKept = Option<(Arc<SequenceBlock>, usize)>;

/// The read of a `-db` stream that failed in the search of one query batch, with the smallest
/// ordinal id (`Subjects::Stream`), or the subject of the batch whose search LOSAT does not
/// support (two HSPs of the one-hit word finder whose chunk merge NCBI leaves undefined,
/// `OneHitChunkMerge::refuse` in `blast_engine/run_impl.rs`). The subjects after it are not
/// searched (`skips`); those before it are, so the error kept does not depend on the threads or
/// their order. Each query batch has its own: a batch that the pipe runs early does not stop the
/// batch before it.
pub(crate) struct XReadFailure {
    /// The ordinal id of the failed read kept, or `usize::MAX`.
    first: AtomicUsize,
    failed: Mutex<Option<(usize, anyhow::Error)>>,
}

impl XReadFailure {
    pub(crate) const fn new() -> Self {
        Self {
            first: AtomicUsize::new(usize::MAX),
            failed: Mutex::new(None),
        }
    }

    /// Whether subject `oid` comes after a read that failed.
    pub(crate) fn skips(&self, oid: usize) -> bool {
        oid > self.first.load(Ordering::Relaxed)
    }

    /// The read of subject `oid` failed with `error`.
    pub(crate) fn fail(&self, oid: usize, error: anyhow::Error) {
        let mut failed = self.failed.lock().unwrap_or_else(|e| e.into_inner());
        if failed.as_ref().is_none_or(|(first, _)| oid < *first) {
            *failed = Some((oid, error));
            self.first.fetch_min(oid, Ordering::Relaxed);
        }
    }

    /// Whether a read failed (without the lock of `take`).
    pub(crate) fn failed(&self) -> bool {
        self.first.load(Ordering::Relaxed) != usize::MAX
    }

    /// The error of the failed read kept, if any, taken.
    pub(crate) fn take(&self) -> Option<anyhow::Error> {
        let mut failed = self.failed.lock().unwrap_or_else(|e| e.into_inner());
        failed.take().map(|(_, error)| error)
    }
}

// No NCBI counterpart: the variable is read once per process; it chooses how LOSAT resolves a
// record's letters, not a value NCBI computes.
/// LOSAT_X_TBXRESOLVE: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
pub(crate) fn x_resolve_by_table() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXRESOLVE"))
}

// No NCBI counterpart: the variable is read once per process; it chooses where LOSAT holds the
// letters of `-subject`, not a value NCBI computes.
/// LOSAT_X_TBXSTORE: a strict-set switch (on unless it is `0`; `utils/xswitch.rs`).
pub(crate) fn x_tbx_store() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXSTORE"))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:424-425,441-450
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
    let mut random = crate::utils::random_gen::NcbiRandom::new(letters.len() as u32);
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:94-103
/// The ncbi4na code of a letter is that of `GeneticCode::get` (both cases, `U` as `T`);
/// another byte is a gap (0).
#[cfg(test)]
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
#[cfg(test)]
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
#[cfg(test)]
const KEPT_BUFFER: usize = 8 << 20;

#[cfg(test)]
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:909-911,974-976
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407,1451-1455
// A subject set (seqsrc_multiseq.cpp:293-298) and a database (seqsrc_seqdb.cpp:183-187) give
// their real totals; `-dbsize` replaces only the effective lengths' database length.
/// The numbers of the subject set that the statistics read, once per run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SubjectSetNumbers {
    /// The bases of every subject (`-db`: the index's total).
    pub(crate) total: i64,
    /// The number of subjects (`-db`: the index's count).
    pub(crate) count: i64,
    /// `BlastSeqSrcGetTotLen`, which the linking cutoffs read (blast_engine.c:1407, 1451-1455):
    /// `total`, except for a database whose alias files give `STATS_TOTLEN` (the statistics' total)
    /// or whose list kept no sequence (`GetTotalLength`, `SeqDb::total_length`).
    pub(crate) link_total: i64,
    /// A database's `BlastSeqSrcGetAvgSeqLen` (`DbStatistics::average`), the subject length of
    /// the word cutoffs (blast_engine.c:1372-1374); `None` for a FASTA set, whose average is
    /// `total / count`.
    pub(crate) average: Option<i32>,
    /// `-dbsize`.
    pub(crate) dbsize: Option<i64>,
    /// `-searchsp`.
    pub(crate) searchsp: Option<i64>,
    /// A database's numbers (`SeqDb::statistics`): the preliminary stage's option spaces and the
    /// statistics of each subject read them (`SubjectSetNumbers::preliminary_spaces`,
    /// `per_subject`); `None` for a FASTA set.
    pub(crate) statistics: Option<crate::seqdb::DbStatistics>,
}

impl SubjectSetNumbers {
    /// The numbers of the records of `-subject`.
    pub(crate) fn of_records(records: &[FastaRecord], dbsize: Option<i64>) -> Self {
        Self {
            total: records.iter().map(|record| record.seq().len() as i64).sum(),
            count: records.len() as i64,
            link_total: records.iter().map(|record| record.seq().len() as i64).sum(),
            average: None,
            dbsize,
            searchsp: None,
            statistics: None,
        }
    }

    /// The numbers with `-searchsp`.
    pub(crate) fn with_searchsp(self, searchsp: Option<i64>) -> Self {
        Self { searchsp, ..self }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:821-846
    /// The length adjustment and the effective search space of a context of `query_len_aa`
    /// residues (`try_compute_eff_lengths_tblastx`), or the option's space other than 0 as the
    /// space (`-searchsp`, or a strategy's `EffectiveSearchSpace`, which can be negative).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:676-696
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1020-1035
    // `stats::setup::eff_lengths_database` and `calc_eff_lengths` (one owner for the six
    // programs, line RG D1a-3): `-dbsize` above 0 or the total, divided by 3 for the translated
    // subject, NCBI's `Int4` count, and a space other than 0 as every context's search space.
    pub(crate) fn context_lengths(
        &self,
        query_len_aa: i64,
        karlin: &crate::stats::KarlinParams,
    ) -> Result<super::ncbi_cutoffs::EffLengthsResult, SearchSpaceOverflow> {
        let searchsp = self.searchsp.unwrap_or(0);
        self.context_lengths_with(query_len_aa, karlin, searchsp)
    }

    /// `context_lengths` with the context's option space `searchsp` (`-searchsp`, or the
    /// preliminary stage's: `preliminary_spaces`).
    pub(crate) fn context_lengths_with(
        &self,
        query_len_aa: i64,
        karlin: &crate::stats::KarlinParams,
        searchsp: i64,
    ) -> Result<super::ncbi_cutoffs::EffLengthsResult, SearchSpaceOverflow> {
        let none = super::ncbi_cutoffs::EffLengthsResult {
            length_adjustment: 0,
            eff_searchsp: searchsp,
        };
        let Some((db_length, count)) = crate::stats::setup::eff_lengths_database(
            self.total,
            self.count as i32,
            self.dbsize.unwrap_or(0),
            searchsp != 0,
            true,
        ) else {
            return Ok(none);
        };
        if query_len_aa <= 0 {
            return Ok(none);
        }
        let lengths = crate::stats::setup::calc_eff_lengths(
            query_len_aa,
            db_length,
            count,
            karlin,
            searchsp,
            "TBLASTX",
        )?;
        Ok(super::ncbi_cutoffs::EffLengthsResult {
            length_adjustment: i64::from(lengths.length_adjustment),
            eff_searchsp: lengths.eff_searchsp,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-220 (cited at
    // `stats::setup::preliminary_option_space`)
    /// The option spaces of the preliminary stage for the contexts `contexts` (each one's query
    /// length in residues, 0 for an invalid one, and its ungapped Karlin block, which the
    /// effective lengths of TBLASTX read): where they differ from the setup's (a database whose
    /// statistics' total before the fallback to OID 0 is not positive) or where the statistics of
    /// each subject read them (`per_subject`); `None` otherwise.
    pub(crate) fn preliminary_spaces<'k>(
        &self,
        contexts: impl Iterator<Item = (i64, &'k crate::stats::KarlinParams)>,
    ) -> Result<Option<Vec<i64>>, SearchSpaceOverflow> {
        let Some(statistics) = self.statistics else {
            return Ok(None);
        };
        if statistics.gumbel_db_length > 0 && !statistics.per_subject() {
            return Ok(None);
        }
        let searchsp = self.searchsp.unwrap_or(0);
        contexts
            .map(|(query_length, karlin)| {
                crate::stats::setup::preliminary_option_space(
                    query_length,
                    statistics.gumbel_db_length,
                    statistics.request_count,
                    self.dbsize.unwrap_or(0),
                    searchsp,
                    true,
                    karlin,
                    "TBLASTX",
                )
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some)
    }

    /// Whether the setup's own spaces are not the preliminary stage's (the statistics' total before
    /// the fallback to OID 0 is not positive: `preliminary_spaces`).
    pub(crate) fn preliminary_spaces_differ(&self) -> bool {
        self.statistics
            .is_some_and(|statistics| statistics.gumbel_db_length <= 0)
    }

    /// Whether NCBI's preliminary search computes the statistics of each subject on its own
    /// (`DbStatistics::per_subject`).
    pub(crate) fn per_subject(&self) -> bool {
        self.statistics
            .is_some_and(|statistics| statistics.per_subject())
    }

    /// The effective lengths of a context of `query_len_aa` residues for the subject of
    /// `subject_nucl` bases (`BLAST_OneSubjectUpdateParameters`, `stats::setup::subject_eff_lengths`),
    /// with the context's preliminary option space `option_space`.
    pub(crate) fn subject_context_lengths(
        &self,
        query_len_aa: i64,
        subject_nucl: i64,
        option_space: i64,
        karlin: &crate::stats::KarlinParams,
    ) -> Result<super::ncbi_cutoffs::EffLengthsResult, SearchSpaceOverflow> {
        let lengths = crate::stats::setup::subject_eff_lengths(
            query_len_aa,
            subject_nucl,
            self.count as i32,
            self.dbsize.unwrap_or(0),
            option_space,
            true,
            karlin,
            "TBLASTX",
        )?;
        Ok(super::ncbi_cutoffs::EffLengthsResult {
            length_adjustment: i64::from(lengths.length_adjustment),
            eff_searchsp: lengths.eff_searchsp,
        })
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
/// LOSAT's refusal of a search whose effective search space overflows NCBI's `Int8`
/// (`db_length` bases of the database, or of `-dbsize`, against a context of `query_length`
/// residues): NCBI's arithmetic overflows there (as BLASTP's and BLASTX's `-dbsize`).
/// LOSAT's refusal of an effective search space past NCBI's `Int8` (the shared one of
/// `stats::setup::calc_eff_lengths`).
pub(crate) use crate::stats::setup::SearchSpaceOverflow;

#[cfg(test)]
#[path = "subjects_tests.rs"]
mod tests;
