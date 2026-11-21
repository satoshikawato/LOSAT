//! The seed search of one MegaBLAST index volume for one query batch (`CDbIndex::Search`):
//! the seeds of the one-hit search (`CSearch<LEGACY, ONE_HIT>`) of every sequence chunk of the
//! volume, kept as NCBI's search results (`CDbIndex::CSearchResults`) for BLASTN's indexed word
//! finder.
//!
//! For each lookup segment of the query (its unmasked intervals, in order) the search walks the
//! Nmers of the segment (`CNmerIterator`), and for each Nmer the stored positions of its list
//! (`CPreOrderedOffsetIterator`: a larger word size visits fewer of them) as roots, kept per
//! local id (`CSeedRoots`). Then, per local id, each root is extended left and right by exact
//! match on the packed sequences and tracked on its diagonal (`CTrackedSeeds<ONE_HIT>`); a seed
//! of at least the word size is saved by the left end of its exact run into the list of its chunk
//! (`SaveSeed`, `MapSubjOff`). Lists keep NCBI's order of saves, which the word finder reads in
//! order. The arithmetic is NCBI's: `TSeqPos` and `TWord` are 32-bit and wrap, `unsigned long`
//! is 64-bit.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:99-1441

use std::sync::atomic::AtomicU32;

use super::volume::IndexVolume;
use super::{code_bits, Word, CR};

/// `TSeqPos`.
type SeqPos = u32;

/// `CDbIndex::SSearchOptions`: the word size and the two-hit window (0 for the one-hit search,
/// the only one here).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SearchOptions {
    pub(crate) word_size: u64,
    pub(crate) two_hits: u64,
}

/// An initial hit that the search saves (`BLAST_SaveInitialHit(hitlist, q_off, s_off, NULL)`):
/// the left end of the seed in the concatenated query and in its chunk.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Seed {
    pub(crate) q_off: u32,
    pub(crate) s_off: u32,
}

/// `CDbIndex::CSearchResults` of a search: the seeds of each chunk of the volume (1-based
/// chunk ids, the sentinel chunk included), and the first chunk of each sequence.
#[derive(Debug)]
pub(crate) struct VolumeSeeds {
    word_size: u64,
    /// `map_`: the first chunk of each sequence (`stop - start` words).
    map: Vec<Word>,
    /// `results_`, by chunk id less one (`start_` is 0).
    results: Vec<Option<Vec<Seed>>>,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex.hpp:497-614
impl VolumeSeeds {
    /// `GetResults(seq)`: the seeds of chunk id `seq`.
    fn by_id(&self, seq: u32) -> Option<&[Seed]> {
        if seq == 0 {
            return None;
        }
        self.results.get(seq as usize - 1)?.as_deref()
    }

    /// `MapSubject(subj, chunk)`: the chunk id of chunk `chunk` of sequence `subj`, 0 for a
    /// sequence past the volume's.
    fn map_subject(&self, subj: u32, chunk: u32) -> u32 {
        match self.map.get(subj as usize) {
            Some(&first) => first.wrapping_add(chunk),
            None => 0,
        }
    }

    /// `GetResults(subj, chunk)`: the seeds of chunk `chunk` of the volume's sequence `subj`
    /// (`None`: none).
    pub(crate) fn seeds(&self, subj: u32, chunk: u32) -> Option<&[Seed]> {
        self.by_id(self.map_subject(subj, chunk))
    }

    /// Whether chunk `chunk` of the volume's sequence `subj` is one of the sequence's chunks in
    /// the index (a sequence past the volume's has none, and NCBI's `GetResults` gives it
    /// nothing). A database sequence with more chunks than the index's sequence (an index not
    /// made from the database) reads the next sequence's seeds, which NCBI's `GetResults` moves
    /// out of the next sequence's list: what that sequence finds depends on the order of the
    /// search threads.
    pub(crate) fn chunk_in_sequence(&self, subj: u32, chunk: u32) -> bool {
        let Some(&first) = self.map.get(subj as usize) else {
            return true;
        };
        let end = match self.map.get(subj as usize + 1) {
            Some(&next) => next,
            None => (self.results.len() as u32).wrapping_add(1),
        };
        u64::from(first) + u64::from(chunk) < u64::from(end)
    }

    /// `CheckResults(subj)`: whether a chunk of the volume's sequence `subj` has seeds; the
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
}

/// An index volume of a search: its path, its first OID and its number of OIDs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IndexedVolume {
    pub(crate) path: String,
    pub(crate) start_oid: u32,
    pub(crate) n_oids: u32,
}

/// The MegaBLAST index that a BLASTN search uses: its volumes in OID order (`CIndexedDb_New`'s
/// `volumes_`, each index volume of each database volume; or `CIndexedDb_Old`'s index volumes,
/// whose OIDs follow one another from 0).
#[derive(Debug)]
pub(crate) struct LoadedIndex {
    pub(crate) old_style: bool,
    pub(crate) volumes: Vec<IndexedVolume>,
    /// The searches of the index so far (`DoPreSearch` calls: one a query batch or query chunk).
    pub(crate) searches: AtomicU32,
}

impl LoadedIndex {
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

/// `SSeedRoot`: a query offset (the last base of its Nmer), a subject offset (with
/// `min_offset` added, or a boundary code before the root that it qualifies), and the query
/// segment.
#[derive(Clone, Copy, Debug)]
struct SeedRoot {
    qoff: SeqPos,
    soff: SeqPos,
    qstart: SeqPos,
    qstop: SeqPos,
}

/// `STrackedSeed<ONE_HIT>`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct TrackedSeed {
    qoff: SeqPos,
    soff: SeqPos,
    len: SeqPos,
    qright: SeqPos,
}

/// `CTrackedSeeds<ONE_HIT>` of one local id: the tracked seeds in diagonal order with the
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

    /// `GetNumChunks(lid)`.
    fn num_chunks(&self, lid: usize) -> usize {
        self.lid_map[4 * lid + 1].wrapping_sub(self.lid_map[4 * lid]) as usize
    }

    /// `SetSubjInfo(lid)`: the local id's first and one past its last byte.
    fn subj_info(&self, lid: usize) -> (Word, Word) {
        (self.lid_map[4 * lid + 2], self.lid_map[4 * lid + 3])
    }

    /// `DecodeOffset(offset)`: the local id and its position (with `min_offset` added).
    fn decode_offset(&self, offset: Word) -> (u32, SeqPos) {
        let offset = offset.wrapping_sub(self.min_offset as Word);
        (
            offset >> self.offset_bits,
            (self.min_offset + u64::from(offset & self.offset_mask) * self.stride) as SeqPos,
        )
    }

    /// `MapSubjOff(lid, soff)`: the chunk of the local id (relative) whose data hold position
    /// `soff` (the last one that starts at or before it) and the position in that chunk.
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

/// `CPreOrderedOffsetIterator`: the offsets of one Nmer's list, pass by pass.
struct OffsetIter {
    /// The index in the offset lists of the word read last (`curr_`), `None` for an empty list.
    curr: Option<usize>,
    offset: Word,
    more: u64,
    modulus: u64,
    boundary: bool,
    min_offset: u64,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/dbindex_sp.hpp:214-276
impl OffsetIter {
    /// `Reset(offset_data, key, ws)`.
    fn new(search: &Search<'_>, key: Word) -> Result<OffsetIter, Corrupt> {
        let h = search.hkey_width - 1;
        let s = search.stride;
        let w = search.ws_hint;
        let mut more = w.wrapping_sub(h) / s;
        let modulus = search.options.word_size.wrapping_sub(h) / s;
        let cache = *search.hash_table.get(key as usize).ok_or(Corrupt)?;
        let curr = if cache != 0 {
            Some(cache as usize - 1)
        } else {
            more = 0;
            None
        };
        Ok(OffsetIter {
            curr,
            offset: 0,
            more,
            modulus,
            boundary: false,
            min_offset: search.min_offset,
        })
    }

    /// `More()`.
    fn more(&self) -> bool {
        self.more != 0
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
            self.more = 0;
            return Ok(false);
        }
        if u64::from(self.offset) < self.min_offset {
            self.boundary = true;
            Ok(true)
        } else if self.boundary {
            self.boundary = false;
            Ok(true)
        } else if u64::from(self.offset) % self.more == 0 {
            Ok(true)
        } else {
            self.more = if self.more <= self.modulus {
                0
            } else {
                self.more - 1
            };
            self.curr = Some(curr - 1);
            Ok(false)
        }
    }
}

/// `CSearch<LEGACY, ONE_HIT>`'s state.
struct Search<'a> {
    hash_table: &'a [Word],
    offsets: &'a [Word],
    seq_store: &'a [u8],
    /// `IndexVolume::seq_store_padding`.
    store_padding: usize,
    map: SubjectMap<'a>,
    hkey_width: u64,
    stride: u64,
    ws_hint: u64,
    min_offset: u64,
    code_bits: u64,
    query: &'a [u8],
    options: SearchOptions,
    seeds: Vec<TrackedSeeds>,
    /// `subject_`, `subj_start_`, `subj_end_`.
    subject: usize,
    subj_start: Word,
    subj_end: Word,
    qoff: SeqPos,
    soff: SeqPos,
    qstart: SeqPos,
    qstop: SeqPos,
    /// `roots_`: the roots of each local id (and of one more), in the order added, and their
    /// total.
    roots: Vec<Vec<SeedRoot>>,
    total_roots: u64,
}

/// `CSeedRoots::LIMIT_ROOTS`: the roots kept before the seeds are computed from them.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:281-283
const LIMIT_ROOTS: u64 = 16 * 1024 * 1024;

// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1411-1437
/// `CDbIndex::Search(query, locs, options)` of the one-hit search: the seeds of `volume` for the
/// query batch `query` (the concatenated query in BLASTNA, NCBI's `query->sequence`) and its
/// lookup segments `segments` (`[start, stop)` in the concatenated query, in NCBI's order).
pub(crate) fn search(
    volume: &IndexVolume,
    query: &[u8],
    segments: &[(u32, u32)],
    options: SearchOptions,
) -> Result<VolumeSeeds, Corrupt> {
    debug_assert_eq!(options.two_hits, 0, "the two-hit search is not ported");
    let map = SubjectMap {
        lid_map: volume.lid_map(),
        chunks: volume.chunks(),
        offset_bits: volume.offset_bits,
        offset_mask: (1u32 << volume.offset_bits) - 1,
        stride: volume.stride,
        min_offset: volume.min_offset(),
    };
    let num_lids = map.num_subjects() - 1;
    let seeds = (0..num_lids)
        .map(|lid| TrackedSeeds {
            seeds: Vec::new(),
            it: 0,
            hitlists: vec![None; map.num_chunks(lid)],
        })
        .collect();
    let mut search = Search {
        hash_table: volume.hash_table(),
        offsets: volume.offsets(),
        seq_store: volume.seq_store(),
        store_padding: volume.seq_store_padding(),
        hkey_width: volume.hkey_width,
        stride: volume.stride,
        ws_hint: volume.ws_hint,
        min_offset: volume.min_offset(),
        code_bits: code_bits(volume.stride),
        roots: vec![Vec::new(); map.num_subjects()],
        map,
        query,
        options,
        seeds,
        subject: 0,
        subj_start: 0,
        subj_end: 0,
        qoff: 0,
        soff: 0,
        qstart: 0,
        qstop: 0,
        total_roots: 0,
    };
    // `operator()`: each segment, then the seeds of the roots left.
    for &(start, stop) in segments {
        search.qstart = start;
        search.qstop = stop;
        search.search_int()?;
    }
    search.compute_seeds()?;
    let chunks = volume.chunks().len();
    let mut results: Vec<Option<Vec<Seed>>> = vec![None; chunks];
    let mut k = 1usize;
    for lid in 0..num_lids {
        search.finalize(lid)?;
        for hitlist in std::mem::take(&mut search.seeds[lid].hitlists) {
            // `SetResults(k, ...)` keeps only the ids of the volume's chunks.
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

impl Search<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:99-180
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1260-1307
    /// `SearchInt()`: the roots of the Nmers of the segment `[qstart, qstop)`; the seeds are
    /// computed (and the roots dropped) whenever more than `LIMIT_ROOTS` are kept, after a
    /// query position.
    fn search_int(&mut self) -> Result<(), Corrupt> {
        let hkey_width = self.hkey_width;
        let hkey_mask: Word = ((1u64 << (2 * hkey_width)) - 1) as Word;
        let (start, stop) = (self.qstart, self.qstop);
        let mut nmer: Word = 0;
        let mut count: u64 = 0;
        let mut pos = start;
        while pos < stop {
            let letter = Word::from(*self.query.get(pos as usize).ok_or(Corrupt)?);
            pos += 1;
            if letter >= 4 {
                count = 0;
                nmer = 0;
                continue;
            }
            nmer = ((nmer << 2) & hkey_mask) + letter;
            count += 1;
            if count < hkey_width {
                continue;
            }
            let mut off_it = OffsetIter::new(self, nmer)?;
            self.qoff = pos - 1;
            let offsets = self.offsets;
            while off_it.more() {
                self.subject = 0;
                while off_it.next(offsets)? {
                    let offset = off_it.offset;
                    if u64::from(offset) < self.min_offset {
                        // The offset that a boundary code qualifies; a list that ends after the
                        // code is corrupt (NCBI's iterator then reads the next list with a pass
                        // of 0, a division by zero).
                        if !off_it.next(offsets)? {
                            return Err(Corrupt);
                        }
                        let soff = self.decode_offset(off_it.offset)?;
                        let r1 = SeedRoot {
                            qoff: self.qoff,
                            soff: offset,
                            qstart: self.qstart,
                            qstop: self.qstop,
                        };
                        let r2 = SeedRoot { soff, ..r1 };
                        let roots = &mut self.roots[self.subject];
                        roots.push(r1);
                        roots.push(r2);
                        self.total_roots += 2;
                    } else {
                        let soff = self.decode_offset(offset)?;
                        self.roots[self.subject].push(SeedRoot {
                            qoff: self.qoff,
                            soff,
                            qstart: self.qstart,
                            qstop: self.qstop,
                        });
                        self.total_roots += 1;
                    }
                }
            }
            if self.total_roots > LIMIT_ROOTS {
                let (old_qstart, old_qstop) = (self.qstart, self.qstop);
                self.compute_seeds()?;
                for roots in &mut self.roots {
                    *roots = Vec::new();
                }
                self.total_roots = 0;
                self.qstart = old_qstart;
                self.qstop = old_qstop;
            }
        }
        Ok(())
    }

    /// `DecodeOffset(offset)` of `CSearch`: the position, with the local id kept as the current
    /// subject (whose roots can be kept: the local ids and one more).
    fn decode_offset(&mut self, offset: Word) -> Result<SeqPos, Corrupt> {
        let (lid, soff) = self.map.decode_offset(offset);
        if lid as usize >= self.roots.len() {
            return Err(Corrupt);
        }
        self.subject = lid as usize;
        Ok(soff)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1224-1257
    /// `ComputeSeeds()`: each local id's roots in order, from the query offset 0.
    fn compute_seeds(&mut self) -> Result<(), Corrupt> {
        for lid in 0..self.seeds.len() {
            self.subject = lid;
            (self.subj_start, self.subj_end) = self.map.subj_info(lid);
            let roots = std::mem::take(&mut self.roots[lid]);
            if !roots.is_empty() {
                self.qoff = 0;
                let mut j = 0;
                while j < roots.len() {
                    j += self.process_root(&roots[j..])?;
                }
            }
            self.seeds[lid].it = 0;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1194-1221
    /// `ProcessRoot(seeds, root)`: a new query offset, or a subject offset below the last one,
    /// starts the tracked seeds over; a boundary code takes the root after it. Returns the
    /// roots read.
    fn process_root(&mut self, roots: &[SeedRoot]) -> Result<usize, Corrupt> {
        let root = roots[0];
        if self.qoff != root.qoff {
            self.seeds[self.subject].it = 0;
            self.qoff = root.qoff;
        } else if u64::from(root.soff) >= self.min_offset && root.soff < self.soff {
            self.seeds[self.subject].it = 0;
        }
        self.qstart = root.qstart;
        self.qstop = root.qstop;
        let min_offset = self.min_offset as SeqPos;
        if u64::from(root.soff) < self.min_offset {
            let real = *roots.get(1).ok_or(Corrupt)?;
            self.process_boundary_offset(real.soff.wrapping_sub(min_offset), root.soff)?;
            self.soff = real.soff;
            Ok(2)
        } else {
            self.process_offset(root.soff.wrapping_sub(min_offset))?;
            self.soff = root.soff;
            Ok(1)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1177-1191
    /// `ProcessOffset(offset)`: a seed of the Nmer at the current query offset and subject
    /// position `offset`, extended both ways when no tracked seed of its diagonal covers it, and
    /// tracked when it reaches the word size.
    fn process_offset(&mut self, offset: Word) -> Result<(), Corrupt> {
        let mut seed = TrackedSeed {
            qoff: self.qoff,
            soff: offset,
            len: self.hkey_width as SeqPos,
            qright: self.qoff,
        };
        if self.eval_and_update(&seed)? {
            self.extend_left(&mut seed, !0)?;
            self.extend_right(&mut seed, !0)?;
            if u64::from(seed.len) >= self.options.word_size {
                self.append_simple(seed);
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:1142-1174
    /// `ProcessBoundaryOffset(offset, bounds)`: a seed next to the end of a valid interval of
    /// the subject, extended at most to the distances that `bounds` codes (left in the high
    /// bits, right in the low `code_bits`; 0 for more than the code holds).
    fn process_boundary_offset(&mut self, offset: Word, bounds: Word) -> Result<(), Corrupt> {
        let nmaxleft = (u64::from(bounds) >> self.code_bits) as SeqPos;
        let nmaxright = (bounds & ((1u32 << self.code_bits) - 1)) as SeqPos;
        let mut seed = TrackedSeed {
            qoff: self.qoff,
            soff: offset,
            len: self.hkey_width as SeqPos,
            qright: self.qoff,
        };
        self.eval_and_update(&seed)?;
        if nmaxleft > 0 {
            self.extend_left(&mut seed, nmaxleft - 1)?;
        } else {
            self.extend_left(&mut seed, !0)?;
        }
        if nmaxright > 0 {
            self.extend_right(&mut seed, nmaxright - 1)?;
        } else {
            self.extend_right(&mut seed, !0)?;
        }
        if nmaxleft > 0 && nmaxright == 0 && u64::from(seed.len) < self.options.word_size {
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
        match self.seq_store.get(at) {
            Some(&byte) => Ok(byte),
            // The zero bytes of the mapping after the end of the file.
            None if at < self.seq_store.len() + self.store_padding => Ok(0),
            None => Err(Corrupt),
        }
    }

    /// The query letter `at`.
    fn letter(&self, at: i64) -> Result<u8, Corrupt> {
        usize::try_from(at)
            .ok()
            .and_then(|at| self.query.get(at))
            .copied()
            .ok_or(Corrupt)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:998-1066
    /// `ExtendLeft(seed, nmax)`: the exact match to the left of the seed's Nmer, at most
    /// `word_size - hkey_width` bases and `nmax`, bounded by the segment's start and the local
    /// id's first byte: the bases of the Nmer's first byte, then whole bytes, then bases.
    fn extend_left(&self, seed: &mut TrackedSeed, nmax: SeqPos) -> Result<(), Corrupt> {
        let cr = CR as i64;
        let hkey_width = self.hkey_width;
        let sstart = i64::from(self.subj_start);
        let first = u64::from(seed.soff).wrapping_sub(hkey_width - 1);
        let mut spos = sstart + (first / CR) as i64;
        let mut qpos = i64::from(seed.qoff) - (hkey_width as i64 - 1);
        let mut incomplete = (first % CR) as u32;
        let qstart = i64::from(self.qstart);
        let limit = self.options.word_size.wrapping_sub(hkey_width);
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
            spos -= 1;
            let mut qbyte: u8 = 0;
            let mut i = 0u32;
            let mut ambig = false;
            while i < CR as u32 {
                qpos -= 1;
                let letter = self.letter(qpos)?;
                qbyte = qbyte.wrapping_add(letter.wrapping_shl(2 * i));
                if letter > 3 {
                    spos += 1;
                    qpos += i64::from(i) + 1;
                    nmax = i;
                    ambig = true;
                    break;
                }
                i += 1;
            }
            if ambig {
                break;
            }
            if sbyte != qbyte {
                spos += 1;
                qpos += i64::from(i);
                break;
            }
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
    /// `ExtendRight(seed, nmax)`: the exact match to the right of the seed's Nmer, at most
    /// `nmax` bases, bounded by the segment's end and the local id's last byte: the bases of the
    /// Nmer's last byte, then whole bytes, then bases.
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
        while u64::from(nmax) >= CR {
            let sbyte = self.store(spos)?;
            spos += 1;
            let mut qbyte: u8 = 0;
            let mut ambig = false;
            for i in 0..CR as u32 {
                let letter = self.letter(qpos)?;
                if letter > 3 {
                    nmax = i;
                    qpos -= i64::from(i);
                    spos -= 1;
                    ambig = true;
                    break;
                }
                qbyte = (qbyte << 2).wrapping_add(letter);
                qpos += 1;
            }
            if ambig {
                break;
            }
            if sbyte != qbyte {
                spos -= 1;
                qpos -= cr;
                break;
            }
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:675-698
    /// `EvalAndUpdate(seed)` of the one-hit search: the tracked seeds before the seed's
    /// diagonal are passed (those that end before its query offset saved and dropped); false
    /// when a tracked seed lies on its diagonal (the cursor then after it).
    fn eval_and_update(&mut self, seed: &TrackedSeed) -> Result<bool, Corrupt> {
        let lid = self.subject;
        loop {
            let tracked = &mut self.seeds[lid];
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:600-604
    /// `AppendSimple(seed)`: tracked before the cursor.
    fn append_simple(&mut self, seed: TrackedSeed) {
        let tracked = &mut self.seeds[self.subject];
        tracked.seeds.insert(tracked.it, seed);
        tracked.it += 1;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:606-640
    /// `Append(seed, word_size)`: a seed on the diagonal of the tracked seed before the cursor
    /// shortens it to the seed's right end (dropping it below the word size) when it ends
    /// further; any other seed of the word size is tracked before the cursor.
    fn append(&mut self, seed: TrackedSeed) {
        let word_size = self.options.word_size;
        let tracked = &mut self.seeds[self.subject];
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
    /// `SaveSeed(seed)` of the current local id: a seed of positive length at the left end of
    /// its run, into the list of its chunk.
    fn save_seed(&mut self, seed: &TrackedSeed) -> Result<(), Corrupt> {
        if seed.len == 0 {
            return Ok(());
        }
        let qoff = seed.qright.wrapping_sub(seed.len).wrapping_add(1);
        let soff = seed.soff.wrapping_sub(seed.qoff.wrapping_sub(qoff));
        let lid = self.subject;
        let (chunk, chunk_off) = self.map.map_subj_off(lid, soff)?;
        let hitlist = self.seeds[lid]
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_search.cpp:666-673
    /// `Finalize()` of local id `lid`: every tracked seed saved, in order.
    fn finalize(&mut self, lid: usize) -> Result<(), Corrupt> {
        self.subject = lid;
        let seeds = std::mem::take(&mut self.seeds[lid].seeds);
        for seed in &seeds {
            self.save_seed(seed)?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "search_tests.rs"]
mod tests;
