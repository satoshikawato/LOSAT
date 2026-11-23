//! The seed search of one MegaBLAST index volume for one query batch (`CDbIndex::Search`):
//! the seeds of the one-hit search (`CSearch<LEGACY, ONE_HIT>`) or of the two-hit search
//! (`CSearch<LEGACY, TWO_HIT>`, `-window_size` above 0) of every sequence chunk of the volume,
//! kept as NCBI's search results (`CDbIndex::CSearchResults`) for BLASTN's indexed word finder.
//!
//! For each lookup segment of the query (its unmasked intervals, in order) the search walks the
//! Nmers of the segment (`CNmerIterator`), and for each Nmer the stored positions of its list
//! (`CPreOrderedOffsetIterator`: a larger word size visits fewer of them) as roots, kept per
//! local id (`CSeedRoots`). Then, per local id, each root is extended left and right by exact
//! match on the packed sequences and tracked on its diagonal (`CTrackedSeeds`); a seed of at
//! least the word size (the two-hit search: one that follows an earlier hit on its diagonal within
//! the window, or of twice the word size) is saved by the left end of its exact run into the list
//! of its chunk (`SaveSeed`, `MapSubjOff`). Lists keep NCBI's order of saves, which the word
//! finder reads in order. The arithmetic is NCBI's: `TSeqPos` and `TWord` are 32-bit and wrap,
//! `unsigned long` is 64-bit.
//!
//! With `-db_hard_mask` (the Owner's answer C of 2026-11-22 to Q-HARDIDX, approved exception 24
//! of PD-LOSAT-NCBI-DEFECTS) a masked letter of the subject matches no query letter: a root whose
//! Nmer holds one is no root, an extension stops before one, and the two-hit search pairs no two
//! hits of a diagonal with a masked letter between them, so that each seed lies in one unmasked
//! stretch and no masked letter decides a seed (NCBI's search reads the index's letters, masked or
//! not). A local id without masked letters is searched as without the option.
//!
//! NCBI's one thread does all of this in turn. LOSAT collects the roots of consecutive query
//! ranges at once and then computes the seeds of each local id at once, on the search's pool
//! ([`search`]); the seeds and their lists are NCBI's at every thread count. While the roots are
//! collected, each Nmer's word of the hash table and then its list are asked of the memory some
//! Nmers before they are read; the extensions compare four query letters with a byte of the packed
//! sequences at once.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:99-1441

use std::sync::atomic::AtomicU32;

use super::volume::IndexVolume;
use super::{code_bits, Word, CR};
use crate::utils::threading::SearchPool;

/// `TSeqPos`.
type SeqPos = u32;

/// `CDbIndex::SSearchOptions`: the word size and the two-hit window (0 for the one-hit search).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SearchOptions {
    pub(crate) word_size: u64,
    pub(crate) two_hits: u64,
}

/// An initial hit that the search saves (`BLAST_SaveInitialHit` of `hitlist`, `q_off`, `s_off` and
/// a null last argument): the left end of the seed in the concatenated query and in its chunk.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Seed {
    pub(crate) q_off: u32,
    pub(crate) s_off: u32,
}

/// `CDbIndex::CSearchResults` of a search: the seeds of each chunk of the volume (1-based
/// chunk ids, the sentinel chunk included), and the first chunk of each sequence.
#[derive(Debug)]
#[cfg_attr(test, derive(PartialEq))]
pub(crate) struct VolumeSeeds {
    word_size: u64,
    /// `map_`: the first chunk of each sequence (`stop - start` words).
    map: Vec<Word>,
    /// `results_`, by chunk id less one (`start_` is 0).
    results: Vec<Option<Vec<Seed>>>,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:497-614
impl VolumeSeeds {
    /// `GetResults` of `seq`: the seeds of chunk id `seq`.
    fn by_id(&self, seq: u32) -> Option<&[Seed]> {
        if seq == 0 {
            return None;
        }
        self.results.get(seq as usize - 1)?.as_deref()
    }

    /// `MapSubject` of `subj` and `chunk`: the chunk id of chunk `chunk` of sequence `subj`, 0 for
    /// a sequence past the volume's.
    fn map_subject(&self, subj: u32, chunk: u32) -> u32 {
        match self.map.get(subj as usize) {
            Some(&first) => first.wrapping_add(chunk),
            None => 0,
        }
    }

    /// `GetResults` of `subj` and `chunk`: the seeds of chunk `chunk` of the volume's sequence
    /// `subj` (`None`: none).
    pub(crate) fn seeds(&self, subj: u32, chunk: u32) -> Option<&[Seed]> {
        self.by_id(self.map_subject(subj, chunk))
    }

    /// Whether chunk `chunk` of the volume's sequence `subj` is a later sequence's chunk in the
    /// index: a database sequence with more chunks than the index's (an index not made from the
    /// database) then reads a later sequence's seeds, which NCBI's `GetResults` moves out of that
    /// sequence's list, so what that sequence finds depends on the order of the search threads.
    /// (Past the last sequence's chunks NCBI's `GetResults` finds no list.)
    pub(crate) fn chunk_of_a_later_sequence(&self, subj: u32, chunk: u32) -> bool {
        match (self.map.get(subj as usize), self.map.get(subj as usize + 1)) {
            (Some(&first), Some(&next)) => u64::from(first) + u64::from(chunk) >= u64::from(next),
            _ => false,
        }
    }

    /// `CheckResults` of `subj`: whether a chunk of the volume's sequence `subj` has seeds; the
    /// chunk ids from the sequence's first to the next sequence's (to the last chunk for the
    /// last sequence, and from id 0 for a sequence without chunks).
    pub(crate) fn has_seeds(&self, subj: u32) -> bool {
        if subj as usize >= self.map.len() {
            return false;
        }
        let start = self.map_subject(subj, 0);
        let mut end = self.map_subject(subj.wrapping_add(1), 0);
        if end == 0 {
            end = (self.results.len() as u32).wrapping_add(1);
        }
        (start..end).any(|seq| self.by_id(seq).is_some())
    }

    /// `GetWordSize()`: the word size of the search.
    pub(crate) fn word_size(&self) -> u64 {
        self.word_size
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:248-259
    /// The results after BLASTN's word finder has moved each list out (`BlastInitHitListMove`
    /// leaves each list empty but in place): the same chunks have lists, of no seeds.
    pub(crate) fn consumed(&self) -> VolumeSeeds {
        VolumeSeeds {
            word_size: self.word_size,
            map: self.map.clone(),
            results: self
                .results
                .iter()
                .map(|list| list.as_ref().map(|_| Vec::new()))
                .collect(),
        }
    }

    /// The results of a volume of `sequences` sequences of one chunk each and no seeds, for the
    /// tests of the search's users.
    #[cfg(test)]
    pub(crate) fn empty_for_tests(sequences: u32) -> VolumeSeeds {
        VolumeSeeds {
            word_size: 28,
            map: (1..=sequences).collect(),
            results: vec![None; sequences as usize],
        }
    }
}

/// A volume of a search (`SVolumeDescriptor`): its path, its first OID, its number of OIDs, and
/// whether it is an index volume (else a database volume without an index, which the search
/// scans).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IndexedVolume {
    pub(crate) path: String,
    pub(crate) start_oid: u32,
    pub(crate) n_oids: u32,
    pub(crate) has_index: bool,
}

/// `SVolResults` of an index volume (`CIndexedDb_New::results_holder_`): how many search threads
/// still need its results (`ref_count`), and the results that a search left for later ones (in
/// the state the word finder left them: [`VolumeSeeds::consumed`]).
#[derive(Debug, Default)]
pub(crate) struct HeldResults {
    pub(crate) ref_count: i64,
    pub(crate) results: Option<std::sync::Arc<VolumeSeeds>>,
}

/// The MegaBLAST index that a BLASTN search uses: its volumes in OID order (`CIndexedDb_New`'s
/// `volumes_`, each index volume of each database volume with an index and each database volume
/// without one; or `CIndexedDb_Old`'s index volumes, whose OIDs follow one another from 0).
#[derive(Debug)]
pub(crate) struct LoadedIndex {
    pub(crate) old_style: bool,
    pub(crate) volumes: Vec<IndexedVolume>,
    /// The searches of the index so far (`DoPreSearch` calls: one a query batch or query chunk
    /// with a valid context).
    pub(crate) searches: AtomicU32,
    /// `results_holder_` of a new-style index, one for each volume, across its searches.
    pub(crate) held: std::sync::Mutex<Vec<HeldResults>>,
}

impl LoadedIndex {
    /// The index of `volumes` (`old_style`: `CIndexedDb_Old`'s; else `CIndexedDb_New`'s), before
    /// its first search.
    pub(crate) fn new(old_style: bool, volumes: Vec<IndexedVolume>) -> Self {
        LoadedIndex {
            old_style,
            held: std::sync::Mutex::new(volumes.iter().map(|_| HeldResults::default()).collect()),
            volumes,
            searches: AtomicU32::new(0),
        }
    }

    /// Whether some database volume has no index (NCBI's `partial` of a new-style index: the
    /// lookup table `eMixedMBLookupTable`, a mixed search).
    pub(crate) fn mixed(&self) -> bool {
        self.volumes.iter().any(|volume| !volume.has_index)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:208-232
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:497-505
    /// The volume of OID `oid` and the OID in it: the last volume that starts at or before it
    /// (`FindVolume`), or for an old-style index the first volume that ends after it
    /// (`LocateIndex`), the first volume with `oid` itself when none does (its `assert` is not in
    /// NCBI's release build).
    pub(crate) fn locate(&self, oid: u32) -> (usize, u32) {
        if self.old_style {
            match self
                .volumes
                .iter()
                .position(|volume| volume.start_oid.wrapping_add(volume.n_oids) > oid)
            {
                Some(at) => (at, oid.wrapping_sub(self.volumes[at].start_oid)),
                None => (0, oid),
            }
        } else {
            let at = self
                .volumes
                .partition_point(|volume| volume.start_oid <= oid)
                .saturating_sub(1);
            (at, oid.wrapping_sub(self.volumes[at].start_oid))
        }
    }
}

/// A corrupt index volume: a list or an offset that points outside the volume's sections (NCBI
/// reads past them).
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct Corrupt;

/// `SSeedRoot`: a query offset (the last base of its Nmer) and a subject offset (with
/// `min_offset` added, or a boundary code before the root that it qualifies). NCBI's root also
/// holds its lookup segment (`qstart`, `qstop`); LOSAT finds the segment again from the query
/// offset when it processes the root ([`LidSearch::process_root`]), the segments being disjoint
/// and in query order, so a root takes 8 bytes instead of 16.
#[derive(Clone, Copy, Debug)]
struct SeedRoot {
    qoff: SeqPos,
    soff: SeqPos,
}

/// `STrackedSeed<ONE_HIT>`, and `STrackedSeed<TWO_HIT>`'s `second_hit_` (the query end of the
/// earlier hit on the seed's diagonal; 0 for none, and in the one-hit search).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct TrackedSeed {
    qoff: SeqPos,
    soff: SeqPos,
    len: SeqPos,
    qright: SeqPos,
    second_hit: SeqPos,
}

impl TrackedSeed {
    /// The `STrackedSeed` constructor of `qoff`, `soff`, `len` and `qright`.
    fn new(qoff: SeqPos, soff: SeqPos, len: SeqPos, qright: SeqPos) -> Self {
        TrackedSeed {
            qoff,
            soff,
            len,
            qright,
            second_hit: 0,
        }
    }
}

/// `CTrackedSeeds<NHITS>` of one local id: the tracked seeds in diagonal order with the
/// cursor `it_` (an index into `seeds`; `std::list` inserts before it and erases at it), and the
/// saved seeds of each chunk of the local id.
struct TrackedSeeds {
    seeds: Vec<TrackedSeed>,
    it: usize,
    hitlists: Vec<Option<Vec<Seed>>>,
}

/// The subject map of a volume as the search reads it (`CSubjectMap`).
struct SubjectMap<'a> {
    lid_map: &'a [Word],
    chunks: &'a [Word],
    offset_bits: u32,
    offset_mask: Word,
    stride: u64,
    min_offset: u64,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:1105-1222
impl SubjectMap<'_> {
    /// `NumSubjects()`: one more than the number of local ids.
    fn num_subjects(&self) -> usize {
        1 + self.lid_map.len() / 4
    }

    /// `GetNumChunks` of `lid`.
    fn num_chunks(&self, lid: usize) -> usize {
        self.lid_map[4 * lid + 1].wrapping_sub(self.lid_map[4 * lid]) as usize
    }

    /// `SetSubjInfo` of `lid`: the local id's first and one past its last byte.
    fn subj_info(&self, lid: usize) -> (Word, Word) {
        (self.lid_map[4 * lid + 2], self.lid_map[4 * lid + 3])
    }

    /// `DecodeOffset` of `offset`: the local id and its position (with `min_offset` added).
    fn decode_offset(&self, offset: Word) -> (u32, SeqPos) {
        let offset = offset.wrapping_sub(self.min_offset as Word);
        (
            offset >> self.offset_bits,
            (self.min_offset + u64::from(offset & self.offset_mask) * self.stride) as SeqPos,
        )
    }

    /// `MapSubjOff` of `lid` and `soff`: the chunk of the local id (relative) whose data hold
    /// position `soff` (the last one that starts at or before it) and the position in that chunk.
    fn map_subj_off(&self, lid: usize, soff: SeqPos) -> Result<(usize, SeqPos), Corrupt> {
        let start = self.lid_map[4 * lid] as usize;
        let end = self.lid_map[4 * lid + 1] as usize;
        let lid_start = self.lid_map[4 * lid + 2];
        let abs_offset = lid_start.wrapping_add(soff / CR as u32);
        let chunks = &self.chunks[start..end];
        let at = chunks.partition_point(|&chunk| chunk <= abs_offset);
        if at == 0 {
            return Err(Corrupt);
        }
        let res = at - 1;
        let chunk_off = u64::from(soff)
            .wrapping_sub(u64::from(chunks[res].wrapping_sub(lid_start)) * CR)
            as SeqPos;
        Ok((res, chunk_off))
    }
}

/// A pass of NCBI's offset iterator (`more_`, an `unsigned long`) with `c` = ⌈2^64 / `d`⌉ for a
/// `d` from 2 to 2^32 - 1 (0 otherwise): whether `d` divides a 32-bit offset (`n % d == 0`)
/// without a division, as `n * c mod 2^64 < c` (Lemire, Kaser and Kurz, "Faster remainder by
/// direct computation", 2019); 1 divides every offset, and a larger `d` only 0. A `d` of 0 ends
/// the list's passes.
#[derive(Clone, Copy, Debug)]
struct Pass {
    d: u64,
    c: u64,
}

impl Pass {
    fn new(d: u64) -> Self {
        let c = if (2..=u64::from(u32::MAX)).contains(&d) {
            u64::MAX / d + 1
        } else {
            0
        };
        Pass { d, c }
    }

    /// Whether the pass divides `n` (`n % d == 0`; `d` above 0).
    fn divides(self, n: Word) -> bool {
        if self.c != 0 {
            u64::from(n).wrapping_mul(self.c) < self.c
        } else {
            self.d == 1 || n == 0
        }
    }
}

/// `CPreOrderedOffsetIterator`: the offsets of one Nmer's list, pass by pass.
struct OffsetIter {
    /// The index in the offset lists of the word read last (`curr_`), `None` for an empty list.
    curr: Option<usize>,
    offset: Word,
    more: Pass,
    modulus: u64,
    boundary: bool,
    min_offset: u64,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex_sp.hpp:214-276
impl OffsetIter {
    /// `Reset(offset_data, key, ws)` of the Nmer whose hash table word is `cache`: the first
    /// pass and the last of the search are the same for every Nmer ([`Shared::first_pass`]).
    fn new(shared: &Shared<'_>, cache: Word) -> OffsetIter {
        let (curr, more) = if cache != 0 {
            (Some(cache as usize - 1), shared.first_pass)
        } else {
            (None, Pass::new(0))
        };
        OffsetIter {
            curr,
            offset: 0,
            more,
            modulus: shared.last_pass,
            boundary: false,
            min_offset: shared.min_offset,
        }
    }

    /// `More()`.
    fn more(&self) -> bool {
        self.more.d != 0
    }

    /// `Next()`: the next offset of the current pass (`offset`), or the end of the pass.
    fn next(&mut self, offsets: &[Word]) -> Result<bool, Corrupt> {
        let Some(curr) = self.curr else {
            return Ok(false);
        };
        let curr = curr + 1;
        self.curr = Some(curr);
        self.offset = *offsets.get(curr).ok_or(Corrupt)?;
        if self.offset == 0 {
            self.more = Pass::new(0);
            return Ok(false);
        }
        if u64::from(self.offset) < self.min_offset {
            self.boundary = true;
            Ok(true)
        } else if self.boundary {
            self.boundary = false;
            Ok(true)
        } else if self.more.divides(self.offset) {
            Ok(true)
        } else {
            self.more = Pass::new(if self.more.d <= self.modulus {
                0
            } else {
                self.more.d - 1
            });
            self.curr = Some(curr - 1);
            Ok(false)
        }
    }
}

/// What every part of a search reads (`CSearch<LEGACY, NHITS>`'s state but the roots and the
/// state of the local id whose roots it processes): the volume, the query batch and its lookup
/// segments, the options, the first and the last pass of every offset list, and with
/// `-db_hard_mask` the masked letters.
struct Shared<'a> {
    hash_table: &'a [Word],
    offsets: &'a [Word],
    seq_store: &'a [u8],
    /// `IndexVolume::seq_store_padding`.
    store_padding: usize,
    map: SubjectMap<'a>,
    num_lids: usize,
    hkey_width: u64,
    stride: u64,
    min_offset: u64,
    code_bits: u64,
    query: &'a [u8],
    /// The lookup segments (`[start, stop)`, disjoint and in query order).
    segments: &'a [(u32, u32)],
    options: SearchOptions,
    /// `CTrackedSeeds<TWO_HIT>`'s `contig_len_`: twice the word size (`unsigned long`).
    contig_len: u64,
    /// `CPreOrderedOffsetIterator`'s first pass (`(ws_hint - (hkey_width - 1)) / stride`) and
    /// its last (`mod_`: `(word_size - (hkey_width - 1)) / stride`).
    first_pass: Pass,
    last_pass: u64,
    /// `-db_hard_mask`: the masked letters of the volume (`[start, end)` of the packed sequences'
    /// letters, sorted and merged), and each local id's slice of them.
    masked: Vec<(u64, u64)>,
    lid_masks: Vec<(usize, usize)>,
}

impl Shared<'_> {
    /// `DecodeOffset` of `offset`, in `CSearch`: the local id and the position; a local id past
    /// the volume's and the one more that NCBI keeps roots of (and never reads) is corrupt.
    fn decode_offset(&self, offset: Word) -> Result<(usize, SeqPos), Corrupt> {
        let (lid, soff) = self.map.decode_offset(offset);
        if lid as usize > self.num_lids {
            return Err(Corrupt);
        }
        Ok((lid as usize, soff))
    }
}

/// The masked letters `masks` of the sequences of `volume` (`search`) as letters of its packed
/// sequences (from each sequence's first chunk's data), sorted and merged, and for each local id
/// of `map` the part of them that its extensions read. A sequence without chunks has no letters
/// there. Kept out of line: inlined, this setup of a volume spilled registers in the search's
/// per-root loop (+0.5 % instructions of a search without masks, D5b-5 A/B).
#[inline(never)]
fn volume_masked_letters(
    volume: &IndexVolume,
    map: &SubjectMap<'_>,
    masks: &[(u32, Vec<(u32, u32)>)],
) -> (Vec<(u64, u64)>, Vec<(usize, usize)>) {
    let (subjects, chunks) = (volume.subjects(), volume.chunks());
    let mut letters: Vec<(u64, u64)> = Vec::new();
    for (seq, ranges) in masks {
        let seq = *seq as usize;
        let Some(&first) = subjects.get(seq) else {
            continue;
        };
        if first == 0 || subjects.get(seq + 1) == Some(&first) {
            continue;
        }
        let Some(&data) = chunks.get(first as usize - 1) else {
            continue;
        };
        let base = u64::from(data) * CR;
        letters.extend(
            ranges
                .iter()
                .filter(|&&(start, end)| start < end)
                .map(|&(start, end)| (base + u64::from(start), base + u64::from(end))),
        );
    }
    letters.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(letters.len());
    for (start, end) in letters {
        match merged.last_mut() {
            Some(last) if start <= last.1 => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    let lid_masks = (0..map.num_subjects() - 1)
        .map(|lid| {
            // The letters that the local id's extensions read (`subj_end_` bounds them loosely).
            let (start, end) = map.subj_info(lid);
            let (from, to) = (u64::from(start) * CR, u64::from(end) * CR);
            let lo = merged.partition_point(|&(_, stop)| stop <= from);
            (lo, merged.partition_point(|&(begin, _)| begin < to).max(lo))
        })
        .collect();
    (merged, lid_masks)
}

/// `CSeedRoots::LIMIT_ROOTS`: the roots kept before the seeds are computed from them.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:281-283
const LIMIT_ROOTS: u64 = 16 * 1024 * 1024;

/// The query ranges of a search's round for each thread of its pool (at most one per query
/// position of the round): the roots of the ranges are collected at once (each range's within its
/// share of half of `LIMIT_ROOTS`), and then the seeds of each local id from them.
const PARTS_PER_THREAD: usize = 4;

/// The search's threads: the pool's when it runs tasks in parallel, else one.
fn pool_threads(pool: &SearchPool<'_>) -> usize {
    if pool.enabled() {
        pool.threads().max(1)
    } else {
        1
    }
}

/// `work` on each of `items` (with its index), as the tasks of `pool` when it has more than one
/// thread (no item's work reads another's), else in turn; the first error by item.
fn each<T: Send>(
    pool: &SearchPool<'_>,
    items: &mut [T],
    work: impl Fn(usize, &mut T) -> Result<(), Corrupt> + Sync + Send,
) -> Result<(), Corrupt> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool_threads(pool) > 1 {
        use rayon::prelude::*;
        return pool.install(|| {
            items
                .par_iter_mut()
                .enumerate()
                .try_for_each(|(at, item)| work(at, item))
        });
    }
    let _ = pool;
    items
        .iter_mut()
        .enumerate()
        .try_for_each(|(at, item)| work(at, item))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1411-1437
/// `CDbIndex::Search` of `query`, `locs` and `options`: the seeds of `volume` (the one-hit
/// search, or the two-hit search of `options.two_hits` above 0) for the query batch `query` (the
/// concatenated query in BLASTNA, NCBI's `query->sequence`) and its lookup segments `segments`
/// (`[start, stop)` in the concatenated query, disjoint and in NCBI's order, which is the query's),
/// with the masked letters `masks` of `-db_hard_mask` (none without it): each sequence of the
/// volume that has some (numbered from 0 in the volume) and its masks (`[start, end)` letters of
/// the sequence).
///
/// NCBI's one thread collects the roots of the query's Nmers (`SearchInt`), and computes the seeds
/// of each local id from its roots in order (`ComputeSeeds`) when it keeps more than `LIMIT_ROOTS`
/// and at the end; what a local id's seeds are does not depend on when they are computed, nor on
/// another local id. LOSAT does the same in rounds on `pool`: the roots of consecutive query ranges
/// are collected at once, each range's into its own lists, and then the seeds of each local id
/// from its roots of the ranges in order; a range that would pass its share of the round ends
/// the round before the Nmer that passes it, and the next round starts there (`search_rounds`).
/// The seeds are NCBI's at every thread count.
pub(crate) fn search(
    volume: &IndexVolume,
    query: &[u8],
    segments: &[(u32, u32)],
    options: SearchOptions,
    masks: &[(u32, Vec<(u32, u32)>)],
    pool: &SearchPool<'_>,
) -> Result<VolumeSeeds, Corrupt> {
    search_with_limit(volume, query, segments, options, masks, pool, LIMIT_ROOTS)
}

/// [`search`] with rounds of at most about `limit` roots (`LIMIT_ROOTS`; any limit gives the
/// same seeds).
fn search_with_limit(
    volume: &IndexVolume,
    query: &[u8],
    segments: &[(u32, u32)],
    options: SearchOptions,
    masks: &[(u32, Vec<(u32, u32)>)],
    pool: &SearchPool<'_>,
    limit: u64,
) -> Result<VolumeSeeds, Corrupt> {
    // NCBI reads the letters of each segment: one past the query is a corrupt search.
    if segments
        .iter()
        .any(|&(start, stop)| start < stop && stop as usize > query.len())
    {
        return Err(Corrupt);
    }
    debug_assert!(segments.windows(2).all(|pair| pair[0].1 <= pair[1].0));
    let map = SubjectMap {
        lid_map: volume.lid_map(),
        chunks: volume.chunks(),
        offset_bits: volume.offset_bits,
        offset_mask: (1u32 << volume.offset_bits) - 1,
        stride: volume.stride,
        min_offset: volume.min_offset(),
    };
    let num_lids = map.num_subjects() - 1;
    let (masked, lid_masks) = volume_masked_letters(volume, &map, masks);
    let mut lids: Vec<TrackedSeeds> = (0..num_lids)
        .map(|lid| TrackedSeeds {
            seeds: Vec::new(),
            it: 0,
            hitlists: vec![None; map.num_chunks(lid)],
        })
        .collect();
    let h = volume.hkey_width.wrapping_sub(1);
    let shared = Shared {
        hash_table: volume.hash_table(),
        offsets: volume.offsets(),
        seq_store: volume.seq_store(),
        store_padding: volume.seq_store_padding(),
        num_lids,
        hkey_width: volume.hkey_width,
        stride: volume.stride,
        min_offset: volume.min_offset(),
        code_bits: code_bits(volume.stride),
        query,
        segments,
        options,
        contig_len: options.word_size.wrapping_mul(2),
        first_pass: Pass::new(volume.ws_hint.wrapping_sub(h) / volume.stride),
        last_pass: options.word_size.wrapping_sub(h) / volume.stride,
        masked,
        lid_masks,
        map,
    };
    search_rounds(&shared, &mut lids, pool, limit)?;
    // `Finalize` of each local id.
    each(pool, &mut lids, |lid, tracked| {
        LidSearch::new(&shared, lid, tracked).finalize()
    })?;
    let chunks = volume.chunks().len();
    let mut results: Vec<Option<Vec<Seed>>> = vec![None; chunks];
    let mut k = 1usize;
    for tracked in &mut lids {
        for hitlist in std::mem::take(&mut tracked.hitlists) {
            // `SetResults` of `k` keeps only the ids of the volume's chunks.
            if let Some(slot) = results.get_mut(k - 1) {
                *slot = hitlist;
            }
            k += 1;
        }
    }
    let n_seqs = volume.stop.wrapping_sub(volume.start) as usize;
    Ok(VolumeSeeds {
        word_size: options.word_size,
        map: volume.subjects()[..n_seqs.min(volume.subjects().len())].to_vec(),
        results,
    })
}

/// One query range of a round: its query positions `[from, to)`, its roots by local id, the query
/// position after the last one read (`to`, or less where the range stopped at its share; `from`
/// for a range that a range before it ended the round without it), the roots that it counted (those
/// of the local id after the last too), and the memory of its lists after its last round, in roots.
struct Part {
    from: SeqPos,
    to: SeqPos,
    next: SeqPos,
    roots: Vec<Vec<SeedRoot>>,
    kept: u64,
    capacity: u64,
}

/// The rounds of [`search`] with `limit` roots (`LIMIT_ROOTS`): the roots of the query's lookup
/// segments by query range, and after each round the seeds of each local id (`lids`) from its
/// roots of the round's ranges in order. A round's ranges (at most one per query position) share
/// half of `limit`: a range stops before the Nmer that takes it past its share, but for the first
/// Nmer with roots of the round's first range, which it keeps whole (as NCBI keeps the Nmer that
/// passes `LIMIT_ROOTS`). So a round holds at most about half of `limit` roots and that one list,
/// and every round reads at least one Nmer. The first round covers the whole query; each round
/// after it, the query positions that held a quarter of `limit` at the density of the round before
/// (the rest of the query after a round without roots), half of what its ranges may hold, so that
/// a range ends the round only where the roots are denser than that. The lists keep their memory
/// from round to round while it is not above `limit` roots in all.
fn search_rounds(
    shared: &Shared<'_>,
    lids: &mut [TrackedSeeds],
    pool: &SearchPool<'_>,
    limit: u64,
) -> Result<(), Corrupt> {
    let (Some(&(first, _)), Some(&(_, end))) = (shared.segments.first(), shared.segments.last())
    else {
        return Ok(());
    };
    let n_parts = PARTS_PER_THREAD * pool_threads(pool);
    let mut parts: Vec<Part> = (0..n_parts)
        .map(|_| Part {
            from: 0,
            to: 0,
            next: 0,
            roots: (0..shared.num_lids).map(|_| Vec::new()).collect(),
            kept: 0,
            capacity: 0,
        })
        .collect();
    let half = (limit / 2).max(1);
    let mut from = first;
    let mut window = u64::from(end.saturating_sub(first));
    while from < end {
        let to = (u64::from(from) + window).min(u64::from(end)) as SeqPos;
        let span = u64::from(to - from);
        let used = (n_parts as u64).min(span) as usize;
        let budget = (half / used as u64).max(1);
        let round = &mut parts[..used];
        for (k, part) in round.iter_mut().enumerate() {
            part.from = from + (span * k as u64 / used as u64) as SeqPos;
            part.to = from + (span * (k as u64 + 1) / used as u64) as SeqPos;
            part.next = part.from;
            part.kept = 0;
        }
        // A range after one that ended the round is not read (the next round reads it).
        let ended = std::sync::atomic::AtomicUsize::new(used);
        each(pool, round, |k, part| {
            if ended.load(std::sync::atomic::Ordering::Relaxed) < k {
                return Ok(());
            }
            (part.next, part.kept) =
                collect_roots(shared, part.from, part.to, budget, k == 0, &mut part.roots)?;
            if part.next < part.to {
                ended.fetch_min(k, std::sync::atomic::Ordering::Relaxed);
            }
            Ok(())
        })?;
        let kept = match round.iter().position(|part| part.next < part.to) {
            Some(k) => k + 1,
            None => used,
        };
        let next = round[kept - 1].next;
        let round = &round[..kept];
        each(pool, lids, |lid, tracked| {
            let mut search = LidSearch::new(shared, lid, tracked);
            for part in round {
                search.process_roots(&part.roots[lid])?;
            }
            search.tracked.it = 0;
            Ok(())
        })?;
        let roots: u64 = round.iter().map(|part| part.kept).sum();
        window = match roots {
            0 => u64::from(end - next),
            _ => ((half / 2).saturating_mul(u64::from(next - from)) / roots).max(1),
        };
        // The lists of the round emptied; the memory of all the ranges' lists kept for the next
        // round while it is not above `limit` roots.
        for part in &mut parts[..used] {
            part.capacity = part.roots.iter().map(|roots| roots.capacity() as u64).sum();
        }
        let free = parts.iter().map(|part| part.capacity).sum::<u64>() > limit;
        let emptied = if free {
            &mut parts[..]
        } else {
            &mut parts[..used]
        };
        for part in emptied {
            for roots in &mut part.roots {
                if free {
                    *roots = Vec::new();
                } else {
                    roots.clear();
                }
            }
            if free {
                part.capacity = 0;
            }
        }
        from = next;
    }
    Ok(())
}

/// The Nmers whose hash table words are asked of the memory at once, and then read at once
/// (the Nmers with a list are kept, whose lists are then asked of the memory), a round of
/// [`collect_roots`] before their lists are read.
const BLOCK: usize = 64;

/// A hint that the word `at` of `words` will be read (x86_64; nothing elsewhere). The address is
/// not read, so it need not be in `words`.
#[inline(always)]
fn prefetch(words: &[Word], at: usize) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not fault, and reads nothing that the program sees.
    unsafe {
        core::arch::x86_64::_mm_prefetch(
            words.as_ptr().wrapping_add(at).cast::<i8>(),
            core::arch::x86_64::_MM_HINT_T0,
        );
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = (words, at);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:99-180
/// `CNmerIterator` over the lookup segments' query positions `[from, to)`: the Nmers that end
/// there (each read from the letters before it in its segment), in order, a block at a time.
struct Nmers<'s> {
    query: &'s [u8],
    segments: &'s [(u32, u32)],
    hkey_width: u64,
    hkey_mask: Word,
    to: SeqPos,
    /// The segment read (`done` past the last before `to`), its first query position whose Nmer
    /// is given, and its end (within `to`).
    segment: usize,
    done: bool,
    emit: SeqPos,
    end: SeqPos,
    pos: SeqPos,
    nmer: Word,
    count: u64,
}

impl<'s> Nmers<'s> {
    fn new(shared: &Shared<'s>, from: SeqPos, to: SeqPos) -> Self {
        let hkey_width = shared.hkey_width;
        let mut nmers = Nmers {
            query: shared.query,
            segments: shared.segments,
            hkey_width,
            hkey_mask: ((1u64 << (2 * hkey_width)) - 1) as Word,
            to,
            segment: shared.segments.partition_point(|&(_, stop)| stop <= from),
            done: false,
            emit: 0,
            end: 0,
            pos: 0,
            nmer: 0,
            count: 0,
        };
        nmers.start_segment(from);
        nmers
    }

    /// Starts the segment `self.segment` (if it starts before `to`) at its first Nmer at or after
    /// `from`: the letters read from the Nmer's first (or from the segment's start), so that every
    /// Nmer read whole ends at `from` or after.
    fn start_segment(&mut self, from: SeqPos) {
        let (start, stop) = match self.segments.get(self.segment) {
            Some(&(start, stop)) if start < self.to => (start, stop),
            _ => {
                (self.done, self.pos, self.end) = (true, 0, 0);
                return;
            }
        };
        let before = self
            .hkey_width
            .saturating_sub(1)
            .min(u64::from(SeqPos::MAX)) as SeqPos;
        self.emit = start.max(from);
        self.end = stop.min(self.to);
        self.pos = start.max(self.emit.saturating_sub(before));
        self.nmer = 0;
        self.count = 0;
    }

    /// The next Nmers into `keys` (their last letter's query offset and their key), up to the
    /// block: how many (fewer only at the end).
    fn fill(&mut self, keys: &mut [(SeqPos, Word); BLOCK], hash_table: &[Word]) -> usize {
        let mut n = 0;
        loop {
            let (end, emit, mask, width) = (self.end, self.emit, self.hkey_mask, self.hkey_width);
            let (mut pos, mut nmer, mut count) = (self.pos, self.nmer, self.count);
            while pos < end && n < BLOCK {
                let letter = Word::from(self.query[pos as usize]);
                pos += 1;
                if letter >= 4 {
                    count = 0;
                    nmer = 0;
                    continue;
                }
                nmer = ((nmer << 2) & mask) + letter;
                count += 1;
                if count >= width {
                    debug_assert!(pos > emit);
                    prefetch(hash_table, nmer as usize);
                    keys[n] = (pos - 1, nmer);
                    n += 1;
                }
            }
            (self.pos, self.nmer, self.count) = (pos, nmer, count);
            if n == BLOCK || self.done {
                return n;
            }
            // The segment is read: the next one, from its start.
            self.segment += 1;
            self.start_segment(0);
        }
    }
}

/// The roots of one Nmer's list in order (`SearchInt`'s loop on the offsets of an Nmer): each
/// offset of the passes as a root, and a boundary code with the offset that it qualifies as two
/// roots of the offset's local id.
struct ListRoots {
    off_it: OffsetIter,
    qoff: SeqPos,
}

impl ListRoots {
    fn new(shared: &Shared<'_>, qoff: SeqPos, cache: Word) -> Self {
        ListRoots {
            off_it: OffsetIter::new(shared, cache),
            qoff,
        }
    }

    /// The next root of the list (with the second root of a boundary code) and its local id.
    #[inline(always)]
    #[allow(clippy::type_complexity)]
    fn next(
        &mut self,
        shared: &Shared<'_>,
    ) -> Result<Option<(usize, SeedRoot, Option<SeedRoot>)>, Corrupt> {
        loop {
            if !self.off_it.more() {
                return Ok(None);
            }
            if self.off_it.next(shared.offsets)? {
                break;
            }
        }
        let (qoff, offset) = (self.qoff, self.off_it.offset);
        if u64::from(offset) < shared.min_offset {
            // The offset that a boundary code qualifies; a list that ends after the code is
            // corrupt (NCBI's iterator then reads the next list with a pass of 0, a division by
            // zero).
            if !self.off_it.next(shared.offsets)? {
                return Err(Corrupt);
            }
            let (lid, soff) = shared.decode_offset(self.off_it.offset)?;
            Ok(Some((
                lid,
                SeedRoot { qoff, soff: offset },
                Some(SeedRoot { qoff, soff }),
            )))
        } else {
            let (lid, soff) = shared.decode_offset(offset)?;
            Ok(Some((lid, SeedRoot { qoff, soff }, None)))
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1260-1307
/// The roots of the Nmers whose last letter is a query position of `[from, to)` in the lookup
/// segments (`SearchInt` of each segment), by local id into `roots`, in NCBI's order, up to
/// `budget` (the roots of the local id after the last, which NCBI keeps and never reads, are
/// counted and not kept): an Nmer whose roots would pass it is not read, unless it is the range's
/// first with roots and `first` (the round's first range), which keeps it whole. Returns the
/// query position after the last one read (`to`, or the position of the Nmer not read, or the one
/// after the first Nmer kept whole past `budget`) and the roots counted. The Nmers go by blocks:
/// a block's hash table words are asked of the memory, then read (the Nmers without a list dropped
/// and the others' lists asked of the memory) a block later, and their lists read a block after
/// that.
fn collect_roots(
    shared: &Shared<'_>,
    from: SeqPos,
    to: SeqPos,
    budget: u64,
    first: bool,
    roots: &mut [Vec<SeedRoot>],
) -> Result<(SeqPos, u64), Corrupt> {
    let mut nmers = Nmers::new(shared, from, to);
    let mut keys = [(0, 0); BLOCK];
    let mut lists = [(0, 0); BLOCK];
    let (mut n_keys, mut n_lists) = (0, 0);
    let mut kept: u64 = 0;
    loop {
        for &(qoff, cache) in &lists[..n_lists] {
            let before = kept;
            let whole = first && before == 0;
            let mut list = ListRoots::new(shared, qoff, cache);
            while let Some((lid, root, second)) = list.next(shared)? {
                if let Some(roots) = roots.get_mut(lid) {
                    roots.push(root);
                    if let Some(second) = second {
                        roots.push(second);
                    }
                }
                kept += 1 + u64::from(second.is_some());
                if kept > budget && !whole {
                    unkeep(shared, qoff, cache, kept - before, roots)?;
                    return Ok((qoff, before));
                }
            }
            if kept > budget {
                return Ok((qoff + 1, kept));
            }
        }
        n_lists = 0;
        for &(qoff, key) in &keys[..n_keys] {
            let cache = *shared.hash_table.get(key as usize).ok_or(Corrupt)?;
            if cache != 0 {
                prefetch(shared.offsets, cache as usize);
                lists[n_lists] = (qoff, cache);
                n_lists += 1;
            }
        }
        n_keys = nmers.fill(&mut keys, shared.hash_table);
        if n_keys == 0 && n_lists == 0 {
            return Ok((to, kept));
        }
    }
}

/// Takes out of `roots` the first `count` roots of the list of the Nmer at `qoff` (whose hash
/// table word is `cache`), which [`collect_roots`] kept last: each the last of its local id's.
#[cold]
fn unkeep(
    shared: &Shared<'_>,
    qoff: SeqPos,
    cache: Word,
    mut count: u64,
    roots: &mut [Vec<SeedRoot>],
) -> Result<(), Corrupt> {
    let mut list = ListRoots::new(shared, qoff, cache);
    while count > 0 {
        let Some((lid, _, second)) = list.next(shared)? else {
            break;
        };
        let n = 1 + usize::from(second.is_some());
        if let Some(roots) = roots.get_mut(lid) {
            roots.truncate(roots.len() - n);
        }
        count -= n as u64;
    }
    Ok(())
}

/// The query letters `at..at + 4` of `query` as a byte of the packed sequences (the first letter
/// in the high bits), and which of them are not bases (`0xFC` in their byte of the little-endian
/// word, 0 for none): the four letters that NCBI's byte loops of the extensions pack one by one.
fn query_byte(query: &[u8], at: i64) -> Result<(u8, u32), Corrupt> {
    let letters = usize::try_from(at)
        .ok()
        .and_then(|at| query.get(at..))
        .and_then(<[u8]>::first_chunk::<4>)
        .ok_or(Corrupt)?;
    let word = u32::from_le_bytes(*letters);
    // Letter `k` (of 2 bits when a base) moves to bits `30 - 2k` and nothing carries into the high
    // byte: the byte is the four letters in order.
    Ok((
        (word.wrapping_mul(0x4010_0401) >> 24) as u8,
        word & 0xFCFC_FCFC,
    ))
}

/// The query letters `at..at + 16` of `query` as four bytes of the packed sequences in a
/// little-endian word (byte `k` the letters `4k..4k + 4`, as [`query_byte`]), when all are bases:
/// the right extension's step over a long exact match (the left one reads at most the word size
/// less the Nmer).
fn query_bytes16(query: &[u8], at: i64) -> Option<u32> {
    let letters = query
        .get(usize::try_from(at).ok()?..)?
        .first_chunk::<16>()?;
    let (lo, hi) = letters.split_at(8);
    let lo = u64::from_le_bytes(lo.try_into().ok()?);
    let hi = u64::from_le_bytes(hi.try_into().ok()?);
    if (lo | hi) & 0xFCFC_FCFC_FCFC_FCFC != 0 {
        return None;
    }
    let pack = |four: u64| (four as u32).wrapping_mul(0x4010_0401) >> 24;
    Some(pack(lo) | pack(lo >> 32) << 8 | pack(hi) << 16 | pack(hi >> 32) << 24)
}

/// One local id's part of a search (`CSearch` while `ComputeSeeds` reads the local id's roots):
/// the local id and its tracked seeds, its bytes of the packed sequences, its masked letters
/// (`-db_hard_mask`: its slice of them and its first letter), the lookup segment of the roots
/// read, and the state that `ProcessRoot` keeps between roots.
struct LidSearch<'s, 'a> {
    shared: &'s Shared<'a>,
    lid: usize,
    tracked: &'s mut TrackedSeeds,
    subj_start: Word,
    subj_end: Word,
    mask_lo: usize,
    mask_hi: usize,
    mask_base: u64,
    /// The index of the lookup segment of the roots read, `qstart_` and `qstop_`.
    segment: usize,
    qstart: SeqPos,
    qstop: SeqPos,
    qoff: SeqPos,
    soff: SeqPos,
}

impl<'s, 'a> LidSearch<'s, 'a> {
    /// The search of local id `lid` with its tracked seeds `tracked`, from the query offset 0
    /// (`ComputeSeeds`; NCBI's other state of the search between local ids is not read before it
    /// is set, the cursor being at the start).
    fn new(shared: &'s Shared<'a>, lid: usize, tracked: &'s mut TrackedSeeds) -> Self {
        let (subj_start, subj_end) = shared.map.subj_info(lid);
        let (mask_lo, mask_hi) = shared.lid_masks[lid];
        LidSearch {
            shared,
            lid,
            tracked,
            subj_start,
            subj_end,
            mask_lo,
            mask_hi,
            mask_base: u64::from(subj_start) * CR,
            segment: 0,
            qstart: 0,
            qstop: 0,
            qoff: 0,
            soff: 0,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1224-1257
    /// `ComputeSeeds()` of the local id's `roots` in order (those of one query range of a round,
    /// after the earlier ranges'), with the masks' search when the local id has masked letters.
    fn process_roots(&mut self, roots: &[SeedRoot]) -> Result<(), Corrupt> {
        let Some(first) = roots.first() else {
            return Ok(());
        };
        if self.qstop <= first.qoff {
            self.segment = self
                .shared
                .segments
                .partition_point(|&(_, stop)| stop <= first.qoff);
        }
        if self.mask_lo < self.mask_hi {
            self.process_roots_of::<true>(roots)
        } else {
            self.process_roots_of::<false>(roots)
        }
    }

    /// Each root of `roots` in order; `MASKED`: the local id has masked letters
    /// (`-db_hard_mask`).
    fn process_roots_of<const MASKED: bool>(&mut self, roots: &[SeedRoot]) -> Result<(), Corrupt> {
        let mut j = 0;
        while j < roots.len() {
            j += self.process_root::<MASKED>(&roots[j..])?;
        }
        Ok(())
    }

    /// With `-db_hard_mask`: `None` when the Nmer whose last letter is `last` (of the current
    /// local id) holds a masked letter, which matches no query letter; else the letters before its
    /// first letter and after its last letter up to the nearest masked letter (`!0`: none on that
    /// side), which its extensions do not pass.
    fn room(&self, last: SeqPos) -> Option<(SeqPos, SeqPos)> {
        let masks = &self.shared.masked[self.mask_lo..self.mask_hi];
        let last = self.mask_base + u64::from(last);
        let first = last.saturating_sub(self.shared.hkey_width - 1);
        let at = masks.partition_point(|&(_, end)| end <= first);
        let right = match masks.get(at) {
            Some(&(start, _)) if start <= last => return None,
            Some(&(start, _)) => SeqPos::try_from(start - last - 1).unwrap_or(!0),
            None => !0,
        };
        let left = match at.checked_sub(1) {
            Some(before) => SeqPos::try_from(first - masks[before].1).unwrap_or(!0),
            None => !0,
        };
        Some((left, right))
    }

    /// With `-db_hard_mask`: whether a masked letter of the current local id lies after letter
    /// `after` and before letter `before`.
    fn masked_between(&self, after: SeqPos, before: SeqPos) -> bool {
        let masks = &self.shared.masked[self.mask_lo..self.mask_hi];
        let after = self.mask_base + u64::from(after) + 1;
        let before = self.mask_base + u64::from(before);
        let at = masks.partition_point(|&(_, end)| end <= after);
        masks
            .get(at)
            .is_some_and(|&(start, _)| start.max(after) < before)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1194-1221
    /// `ProcessRoot` of `seeds` and `root`: a new query offset, or a subject offset below the last
    /// one, starts the tracked seeds over; a boundary code takes the root after it. Returns the
    /// roots read. The root's lookup segment is the one that holds its query offset. `MASKED`: a
    /// root whose Nmer holds a masked letter is passed over as if the list did not hold it, and the
    /// others' extensions stop before the masked letters.
    fn process_root<const MASKED: bool>(&mut self, roots: &[SeedRoot]) -> Result<usize, Corrupt> {
        let root = roots[0];
        let min_offset = self.shared.min_offset;
        let mut room: (SeqPos, SeqPos) = (!0, !0);
        if MASKED {
            let (read, real) = if u64::from(root.soff) < min_offset {
                (2, *roots.get(1).ok_or(Corrupt)?)
            } else {
                (1, root)
            };
            match self.room(real.soff.wrapping_sub(min_offset as SeqPos)) {
                Some(found) => room = found,
                None => return Ok(read),
            }
        }
        if self.qoff != root.qoff {
            self.tracked.it = 0;
            self.qoff = root.qoff;
        } else if u64::from(root.soff) >= min_offset && root.soff < self.soff {
            self.tracked.it = 0;
        }
        if self.qstop <= root.qoff {
            let segments = self.shared.segments;
            while segments.get(self.segment).ok_or(Corrupt)?.1 <= root.qoff {
                self.segment += 1;
            }
            (self.qstart, self.qstop) = segments[self.segment];
        }
        let min_offset = min_offset as SeqPos;
        if u64::from(root.soff) < self.shared.min_offset {
            let real = *roots.get(1).ok_or(Corrupt)?;
            self.process_boundary_offset::<MASKED>(
                real.soff.wrapping_sub(min_offset),
                root.soff,
                room,
            )?;
            self.soff = real.soff;
            Ok(2)
        } else {
            self.process_offset::<MASKED>(root.soff.wrapping_sub(min_offset), room)?;
            self.soff = root.soff;
            Ok(1)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1177-1191
    /// `ProcessOffset` of `offset`: a seed of the Nmer at the current query offset and subject
    /// position `offset`, extended both ways when no tracked seed of its diagonal covers it, and
    /// tracked when it reaches the word size; the extensions at most `room` (`process_root`).
    fn process_offset<const MASKED: bool>(
        &mut self,
        offset: Word,
        room: (SeqPos, SeqPos),
    ) -> Result<(), Corrupt> {
        let hkey_width = self.shared.hkey_width as SeqPos;
        let mut seed = TrackedSeed::new(self.qoff, offset, hkey_width, self.qoff);
        if self.eval_and_update::<MASKED>(&mut seed)? {
            self.extend_left(&mut seed, room.0)?;
            self.extend_right(&mut seed, room.1)?;
            if u64::from(seed.len) >= self.shared.options.word_size {
                self.append_simple(seed);
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1142-1174
    /// `ProcessBoundaryOffset` of `offset` and `bounds`: a seed next to the end of a valid interval
    /// of the subject, extended at most to the distances that `bounds` codes (left in the high
    /// bits, right in the low `code_bits`; 0 for more than the code holds) and `room`
    /// (`process_root`).
    fn process_boundary_offset<const MASKED: bool>(
        &mut self,
        offset: Word,
        bounds: Word,
        room: (SeqPos, SeqPos),
    ) -> Result<(), Corrupt> {
        let code_bits = self.shared.code_bits;
        let nmaxleft = (u64::from(bounds) >> code_bits) as SeqPos;
        let nmaxright = (bounds & ((1u32 << code_bits) - 1)) as SeqPos;
        let hkey_width = self.shared.hkey_width as SeqPos;
        let mut seed = TrackedSeed::new(self.qoff, offset, hkey_width, self.qoff);
        self.eval_and_update::<MASKED>(&mut seed)?;
        if nmaxleft > 0 {
            self.extend_left(&mut seed, (nmaxleft - 1).min(room.0))?;
        } else {
            self.extend_left(&mut seed, room.0)?;
        }
        if nmaxright > 0 {
            self.extend_right(&mut seed, (nmaxright - 1).min(room.1))?;
        } else {
            self.extend_right(&mut seed, room.1)?;
        }
        if nmaxleft > 0 && nmaxright == 0 && u64::from(seed.len) < self.shared.options.word_size {
            seed.len = 0;
            self.append_simple(seed);
        } else {
            self.append(seed);
        }
        Ok(())
    }

    /// The packed byte `at` of the sequences.
    fn store(&self, at: i64) -> Result<u8, Corrupt> {
        let at = usize::try_from(at).map_err(|_| Corrupt)?;
        match self.shared.seq_store.get(at) {
            Some(&byte) => Ok(byte),
            // The zero bytes of the mapping after the end of the file.
            None if at < self.shared.seq_store.len() + self.shared.store_padding => Ok(0),
            None => Err(Corrupt),
        }
    }

    /// The packed bytes `at..at + 4` of the sequences as a little-endian word, where they all lie
    /// in the sequences (`None` before or past them, which the byte loops read one by one).
    fn store4(&self, at: i64) -> Option<u32> {
        let at = usize::try_from(at).ok()?;
        let bytes = self.shared.seq_store.get(at..)?.first_chunk::<4>()?;
        Some(u32::from_le_bytes(*bytes))
    }

    /// The query letter `at`.
    fn letter(&self, at: i64) -> Result<u8, Corrupt> {
        usize::try_from(at)
            .ok()
            .and_then(|at| self.shared.query.get(at))
            .copied()
            .ok_or(Corrupt)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:998-1066
    /// `ExtendLeft` of `seed` and `nmax`: the exact match to the left of the seed's Nmer, at most
    /// `word_size - hkey_width` bases and `nmax`, bounded by the segment's start and the local
    /// id's first byte: the bases of the Nmer's first byte, then whole bytes (four query letters
    /// compared at once; NCBI packs them one by one and stops at the first letter that is not a
    /// base, whose letters before it the base loop then compares), then bases.
    fn extend_left(&self, seed: &mut TrackedSeed, nmax: SeqPos) -> Result<(), Corrupt> {
        let cr = CR as i64;
        let hkey_width = self.shared.hkey_width;
        let sstart = i64::from(self.subj_start);
        let first = u64::from(seed.soff).wrapping_sub(hkey_width - 1);
        let mut spos = sstart + (first / CR) as i64;
        let mut qpos = i64::from(seed.qoff) - (hkey_width as i64 - 1);
        let mut incomplete = (first % CR) as u32;
        let qstart = i64::from(self.qstart);
        let limit = self.shared.options.word_size.wrapping_sub(hkey_width);
        let mut nmax = if u64::from(nmax) < limit {
            nmax
        } else {
            limit as SeqPos
        };
        while nmax > 0 && incomplete > 0 && qpos > qstart {
            let shift = 2 * (CR as u32 - incomplete);
            incomplete -= 1;
            let sbyte = (self.store(spos)? >> shift) & 0x3;
            qpos -= 1;
            if self.letter(qpos)? != sbyte {
                return Ok(());
            }
            nmax -= 1;
            seed.len = seed.len.wrapping_add(1);
        }
        nmax = nmax.min((qpos - qstart) as SeqPos);
        nmax = nmax.min((cr * (spos - sstart)) as SeqPos);
        spos -= 1;
        while u64::from(nmax) >= CR {
            let sbyte = self.store(spos)?;
            let (qbyte, ambiguous) = query_byte(self.shared.query, qpos - cr)?;
            if ambiguous != 0 {
                // The letters after the last that is not a base (read from the right).
                nmax = ambiguous.leading_zeros() / 8;
                break;
            }
            if sbyte != qbyte {
                break;
            }
            spos -= 1;
            qpos -= cr;
            nmax -= CR as SeqPos;
            seed.len = seed.len.wrapping_add(CR as SeqPos);
        }
        let mut i = 0u32;
        while nmax > 0 {
            let sbyte = (self.store(spos)? >> (2 * i)) & 0x3;
            i += 1;
            qpos -= 1;
            if sbyte != self.letter(qpos)? {
                return Ok(());
            }
            seed.len = seed.len.wrapping_add(1);
            nmax -= 1;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1070-1139
    /// `ExtendRight` of `seed` and `nmax`: the exact match to the right of the seed's Nmer, at most
    /// `nmax` bases, bounded by the segment's end and the local id's last byte: the bases of the
    /// Nmer's last byte, then whole bytes (four bytes at once while all their letters match, then
    /// one, four query letters compared at once as in `extend_left`), then bases.
    fn extend_right(&self, seed: &mut TrackedSeed, mut nmax: SeqPos) -> Result<(), Corrupt> {
        let cr = CR as i64;
        let send = i64::from(self.subj_end);
        let mut spos = i64::from(self.subj_start) + i64::from(seed.soff / CR as u32);
        let qend = i64::from(self.qstop);
        let mut qpos = i64::from(seed.qoff) + 1;
        let mut incomplete = seed.soff % CR as u32;
        while nmax > 0 {
            incomplete += 1;
            if incomplete.is_multiple_of(CR as u32) || qpos >= qend {
                break;
            }
            let sbyte = (self.store(spos)? >> (6 - 2 * incomplete)) & 0x3;
            let letter = self.letter(qpos)?;
            qpos += 1;
            if letter != sbyte {
                return Ok(());
            }
            seed.len = seed.len.wrapping_add(1);
            seed.qright = seed.qright.wrapping_add(1);
            nmax -= 1;
        }
        spos += 1;
        nmax = nmax.min((qend - qpos) as SeqPos);
        // `(send - spos) * CR` is an `unsigned long`: a negative difference bounds nothing.
        let store_bases = ((send - spos) as u64).wrapping_mul(CR);
        if u64::from(nmax) > store_bases {
            nmax = store_bases as SeqPos;
        }
        // Four bytes at once while all their letters match (as four turns of the byte loop).
        while u64::from(nmax) >= 4 * CR {
            match (self.store4(spos), query_bytes16(self.shared.query, qpos)) {
                (Some(bytes), Some(letters)) if bytes == letters => {}
                _ => break,
            }
            spos += 4;
            qpos += 4 * cr;
            seed.len = seed.len.wrapping_add(4 * CR as SeqPos);
            seed.qright = seed.qright.wrapping_add(4 * CR as SeqPos);
            nmax -= 4 * CR as SeqPos;
        }
        while u64::from(nmax) >= CR {
            let sbyte = self.store(spos)?;
            let (qbyte, ambiguous) = query_byte(self.shared.query, qpos)?;
            if ambiguous != 0 {
                // The letters before the first that is not a base.
                nmax = ambiguous.trailing_zeros() / 8;
                break;
            }
            if sbyte != qbyte {
                break;
            }
            spos += 1;
            qpos += cr;
            seed.len = seed.len.wrapping_add(CR as SeqPos);
            seed.qright = seed.qright.wrapping_add(CR as SeqPos);
            nmax -= CR as SeqPos;
        }
        let mut i = 2 * (CR as u32 - 1);
        while nmax > 0 {
            nmax -= 1;
            let sbyte = (self.store(spos)? >> i) & 0x3;
            let letter = self.letter(qpos)?;
            qpos += 1;
            if sbyte != letter {
                break;
            }
            seed.len = seed.len.wrapping_add(1);
            seed.qright = seed.qright.wrapping_add(1);
            i = i.wrapping_sub(2);
        }
        Ok(())
    }

    /// `EvalAndUpdate` of `seed` for the search's kind of tracked seeds.
    fn eval_and_update<const MASKED: bool>(
        &mut self,
        seed: &mut TrackedSeed,
    ) -> Result<bool, Corrupt> {
        if self.shared.options.two_hits == 0 {
            self.eval_and_update_one_hit(seed)
        } else {
            self.eval_and_update_two_hits::<MASKED>(seed)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:703-727
    /// `CTrackedSeeds<ONE_HIT>::EvalAndUpdate` of `seed`: the tracked seeds before the seed's
    /// diagonal are passed (those that end before its query offset saved and dropped); false
    /// when a tracked seed lies on its diagonal (the cursor then after it).
    fn eval_and_update_one_hit(&mut self, seed: &TrackedSeed) -> Result<bool, Corrupt> {
        loop {
            let tracked = &mut *self.tracked;
            let Some(&it) = tracked.seeds.get(tracked.it) else {
                return Ok(true);
            };
            let step = seed.qoff.wrapping_sub(it.qoff);
            let it_soff_corr = it.soff.wrapping_add(step);
            if it_soff_corr > seed.soff {
                return Ok(true);
            }
            if it.qright < seed.qoff {
                let at = tracked.it;
                tracked.seeds.remove(at);
                self.save_seed(&it)?;
            } else {
                tracked.it += 1;
                if it_soff_corr == seed.soff {
                    return Ok(false);
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:815-845
    /// `CTrackedSeeds<TWO_HIT>::EvalAndUpdate` of `seed`: the tracked seeds before the seed's
    /// diagonal are passed; one that ends more than the seed's length, the window and three
    /// strides before the seed's end is checked (saved if it qualifies) and dropped; one that
    /// ends before the seed's query offset is dropped if it qualifies and is saved, else it gives
    /// the seed its earlier hit when it lies on the seed's diagonal with a positive length; false
    /// when a tracked seed that reaches the seed's query offset lies on its diagonal (the cursor
    /// then after it). The sums are NCBI's: `TSeqPos` wraps at 32 bits before the `unsigned
    /// long` window and stride are added. `MASKED`: no earlier hit with a masked letter between
    /// its end and the seed's Nmer is the seed's (`-db_hard_mask`, a hit of another stretch).
    fn eval_and_update_two_hits<const MASKED: bool>(
        &mut self,
        seed: &mut TrackedSeed,
    ) -> Result<bool, Corrupt> {
        let reach = self
            .shared
            .options
            .two_hits
            .wrapping_add(self.shared.stride.wrapping_mul(3));
        let hkey_width = self.shared.hkey_width as SeqPos;
        loop {
            let tracked = &mut *self.tracked;
            let Some(&it) = tracked.seeds.get(tracked.it) else {
                return Ok(true);
            };
            let step = seed.qoff.wrapping_sub(it.qoff);
            let it_soff_corr = it.soff.wrapping_add(step);
            if it_soff_corr > seed.soff {
                return Ok(true);
            }
            if u64::from(it.qright.wrapping_add(seed.len)).wrapping_add(reach)
                < u64::from(seed.qright)
            {
                let at = tracked.it;
                tracked.seeds.remove(at);
                self.check_and_save_seed(&it)?;
            } else if it.qright < seed.qoff {
                if self.check_and_save_seed(&it)? {
                    let tracked = &mut *self.tracked;
                    let at = tracked.it;
                    tracked.seeds.remove(at);
                } else {
                    if it_soff_corr == seed.soff
                        && it.len > 0
                        && !(MASKED
                            && self.masked_between(
                                it.soff.wrapping_add(it.qright.wrapping_sub(it.qoff)),
                                seed.soff.wrapping_sub(hkey_width - 1),
                            ))
                    {
                        seed.second_hit = it.qright;
                    }
                    self.tracked.it += 1;
                }
            } else {
                tracked.it += 1;
                if it_soff_corr == seed.soff {
                    return Ok(false);
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:790-801
    /// `CheckAndSaveSeed` of `seed` in the two-hit search: a seed whose earlier hit ends at least
    /// its length, and at most its length and the window, before its end, or a seed of at least
    /// twice the word size, is saved (true).
    fn check_and_save_seed(&mut self, seed: &TrackedSeed) -> Result<bool, Corrupt> {
        let after_hit = seed.second_hit.wrapping_add(seed.len);
        if (seed.second_hit > 0
            && seed.qright >= after_hit
            && u64::from(seed.qright)
                <= u64::from(after_hit).wrapping_add(self.shared.options.two_hits))
            || u64::from(seed.len) >= self.shared.contig_len
        {
            self.save_seed(seed)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:600-604
    /// `AppendSimple` of `seed`: tracked before the cursor.
    fn append_simple(&mut self, seed: TrackedSeed) {
        let tracked = &mut *self.tracked;
        tracked.seeds.insert(tracked.it, seed);
        tracked.it += 1;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:606-640
    /// `Append` of `seed` and `word_size`: a seed on the diagonal of the tracked seed before the
    /// cursor shortens it to the seed's right end (dropping it below the word size) when it ends
    /// further; any other seed of the word size is tracked before the cursor.
    fn append(&mut self, seed: TrackedSeed) {
        let word_size = self.shared.options.word_size;
        let tracked = &mut *self.tracked;
        if tracked.it > 0 {
            let prev = tracked.it - 1;
            let tmp = &mut tracked.seeds[prev];
            let step = seed.qoff.wrapping_sub(tmp.qoff);
            let bs_soff_corr = tmp.soff.wrapping_add(step);
            if bs_soff_corr == seed.soff {
                if seed.qright < tmp.qright {
                    if tmp.len > 0 {
                        tmp.len = tmp.len.wrapping_sub(tmp.qright.wrapping_sub(seed.qright));
                    }
                    if u64::from(tmp.len) < word_size {
                        tracked.seeds.remove(prev);
                        tracked.it -= 1;
                    } else {
                        tmp.qright = seed.qright;
                    }
                }
            } else if u64::from(seed.len) >= word_size {
                tracked.seeds.insert(tracked.it, seed);
                tracked.it += 1;
            }
        } else if u64::from(seed.len) >= word_size {
            tracked.seeds.insert(tracked.it, seed);
            tracked.it += 1;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:578-600
    /// `SaveSeed` of `seed` for the local id: a seed of positive length at the left end of its
    /// run, into the list of its chunk.
    fn save_seed(&mut self, seed: &TrackedSeed) -> Result<(), Corrupt> {
        if seed.len == 0 {
            return Ok(());
        }
        let qoff = seed.qright.wrapping_sub(seed.len).wrapping_add(1);
        let soff = seed.soff.wrapping_sub(seed.qoff.wrapping_sub(qoff));
        let (chunk, chunk_off) = self.shared.map.map_subj_off(self.lid, soff)?;
        let hitlist = self
            .tracked
            .hitlists
            .get_mut(chunk)
            .ok_or(Corrupt)?
            .get_or_insert_with(Vec::new);
        hitlist.push(Seed {
            q_off: qoff,
            s_off: chunk_off,
        });
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:693-699
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:805-811
    /// `Finalize()` of the local id: every tracked seed in order saved (the one-hit search) or
    /// checked and saved if it qualifies (the two-hit search).
    fn finalize(mut self) -> Result<(), Corrupt> {
        let seeds = std::mem::take(&mut self.tracked.seeds);
        for seed in &seeds {
            if self.shared.options.two_hits == 0 {
                self.save_seed(seed)?;
            } else {
                self.check_and_save_seed(seed)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "search_tests.rs"]
mod tests;
