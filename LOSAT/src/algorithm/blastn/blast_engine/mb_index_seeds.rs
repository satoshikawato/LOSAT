//! The search of a query batch with a MegaBLAST database index (`-use_index true` with an index
//! that NCBI's loader kept, [`super::mb_index`]): the index's seeds of each index volume, searched
//! when the preliminary stage first reaches a sequence of the volume (`CIndexedDb_New::UpdateIndex`;
//! every volume of an old-style index the same way, whose seeds NCBI searches before the
//! subjects, with the same result), the sequences without seeds skipped (`CheckOid`), and the word
//! finder of a chunk with seeds (`MB_IndexedWordFinder`): each seed in the order of its list,
//! extended on its diagonal unless an extension on the diagonal already passed its end, and kept
//! with its extension when it reaches the context's cutoff. The chunk's HSP list is then that of
//! every chunk (`prelim_subject::chunk_hsp_list`), without the count of the saved initial hits:
//! NCBI's indexed word finder does not update the ungapped statistics, so the hits of an indexed
//! volume do not size the next query batch. The sequences of a database volume without an index
//! (a mixed search) are scanned with the lookup table as without an index.
//!
//! With `-db_hard_mask` NCBI's word finder takes the list of the index's chunk whose number is the
//! hard chunk's, of another range or of the next sequence; LOSAT searches the index with the masked
//! letters matching nothing (`dbindex::search`) and places the seeds of each sequence into its hard
//! ranges ([`HardSeeds`], the Owner's answer C to Q-HARDIDX, approved exception 24 of
//! PD-LOSAT-NCBI-DEFECTS).
//!
//! NCBI keeps an index volume's results while a search thread may still need them
//! (`results_holder_`, counted by the threads that have not passed the volume), and searches the
//! volume again only once its count is down to 0. A thread that goes on to a database volume
//! without an index does not count down the volumes it passed (`UpdateIndex` returns before), so
//! in a mixed search the last index volume that a search reaches before each such volume that it
//! reaches (by the database's OID list) keeps its results: the next search (a
//! later query batch or query chunk) takes them as the word finder left them, every list moved
//! out and empty (`BlastInitHitListMove`), and finds no seed there; a sequence whose list the
//! earlier search had is read and searched with none, the others skipped. LOSAT reproduces this
//! with NCBI's counts in one thread's order of the sequences ([`plan_search`]). NCBI's search
//! thread then holds such an empty list as its list of initial hits, of no memory: the next
//! initial hit that the scan of a volume without an index saves takes `realloc` of 0 bytes (glibc
//! gives a block of one hit) and the second frees it, after which NCBI reads and frees freed
//! memory (it crashes, or not): LOSAT refuses such a search ([`IndexSeeds::check`]). With more
//! than one thread, which volumes keep their results depends on which thread passes which volume
//! (NCBI crashes in some runs where its one thread does not): LOSAT's search at every thread count
//! is the one-thread search, whose output is NCBI's `-num_threads 1` output (Owner 2026-10-13,
//! 「DB-3: NCBI の 1 スレッド。LOSAT はどのスレッド数でも同じ出力」).

use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use super::prelim_subject::{chunk_hsp_list, PrelimSubject, ScanBatch};
use super::{
    extend_hit_ungapped_approx_ncbi, x_push_ungapped, GapAlignScratch, PrelimHit, QueryContext,
    SeqRanges, SubjectChunk, SubjectScratch, UngappedHit, COMPRESSION_RATIO, MAX_DBSEQ_LEN,
};
use crate::algorithm::blastn::lookup::build_unmasked_ranges;
use crate::dbindex::search::{
    search, HeldResults, IndexedVolume, LoadedIndex, SearchOptions, Seed, VolumeSeeds,
};
use crate::dbindex::volume::IndexVolume;
use crate::seqdb::OidMap;

/// What a search of the index does with each volume (`CIndexedDb_New::UpdateIndex` over the
/// database sequences in order).
#[derive(Debug)]
pub(super) enum VolumePlan {
    /// A database volume without an index: its sequences are scanned.
    Scanned,
    /// An index volume searched for this search; `kept`: its results stay for later searches.
    Searched { kept: bool },
    /// An index volume whose results of an earlier search stay, as the word finder left them.
    Earlier(Arc<VolumeSeeds>),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:812-857
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:877-888
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1364-1366
/// What one search thread of a new-style index does with each volume when it reaches the
/// database sequences of the `visited` volumes in order (`UpdateIndex` at the first sequence of
/// each, then `EndSearchIndication`), with the counts of `held` from the searches before, which it
/// updates: an index volume whose count is not above 0 is searched (its count up by one), one
/// whose count is is taken as it is; reaching an index volume counts down each volume from the
/// last one reached (the first without one reached) to it, reaching a database volume without an
/// index counts none down, and the end counts down each from the last one reached; a count down
/// to 0 drops its results. A searched volume whose count stays above 0 keeps its results for
/// later searches (`kept`, which [`IndexSeeds::check`] stores).
pub(super) fn plan_search(
    volumes: &[IndexedVolume],
    visited: &[bool],
    held: &mut [HeldResults],
) -> Vec<VolumePlan> {
    let mut plan: Vec<VolumePlan> = volumes
        .iter()
        .map(|volume| {
            if volume.has_index {
                VolumePlan::Searched { kept: false }
            } else {
                VolumePlan::Scanned
            }
        })
        .collect();
    let mut searched = vec![false; volumes.len()];
    let count_down = |held: &mut [HeldResults], at: usize| {
        held[at].ref_count -= 1;
        if held[at].ref_count == 0 {
            held[at].results = None;
        }
    };
    let mut last: Option<usize> = None;
    for at in (0..volumes.len()).filter(|&at| visited[at]) {
        if !volumes[at].has_index {
            last = Some(at);
            continue;
        }
        if held[at].ref_count <= 0 {
            held[at].ref_count += 1;
            held[at].results = None;
            searched[at] = true;
        } else {
            let earlier = held[at]
                .results
                .clone()
                .expect("a volume whose count is above 0 holds the results of a search");
            plan[at] = VolumePlan::Earlier(earlier);
        }
        for passed in last.unwrap_or(0)..at {
            count_down(held, passed);
        }
        last = Some(at);
    }
    for passed in last.unwrap_or(0)..volumes.len() {
        count_down(held, passed);
    }
    for (at, searched) in searched.into_iter().enumerate() {
        if searched {
            plan[at] = VolumePlan::Searched {
                kept: held[at].ref_count > 0,
            };
        }
    }
    plan
}

/// The index's seeds of one query batch (`CIndexedDb_New` after `DoPreSearch`, or
/// `CIndexedDb_Old`): the index, the OIDs of the database's ordinal ids, the batch's
/// concatenated query and lookup segments, the search options, the number of the search among
/// the index's searches (from 1; 0 for a batch without a valid context, which NCBI does not
/// search), what the search does with each volume and whether it reaches the volume (a database
/// sequence of it is searched), each index volume's seeds once searched (or `None` after an
/// error), whether the word finder took a list of each index volume, the initial hits (up to 2)
/// that the scan of each database volume without an index saved, and LOSAT's refusal of the
/// search, the first by volume and kind.
pub(super) struct IndexSeeds<'a> {
    index: &'a LoadedIndex,
    oids: OidMap<'a>,
    query: &'a [u8],
    segments: Vec<(u32, u32)>,
    options: SearchOptions,
    search: u32,
    plan: Vec<VolumePlan>,
    visited: Vec<bool>,
    volumes: Vec<OnceLock<Option<Arc<VolumeSeeds>>>>,
    moved: Vec<AtomicBool>,
    scanned_hits: Vec<AtomicU32>,
    error: Mutex<Option<((usize, u8), String)>>,
    /// `-db_hard_mask`: the database and the algorithm whose masked letters the search of each
    /// volume reads.
    hard_masks: Option<(&'a super::database::BlastnDatabase, i32)>,
}

impl<'a> IndexSeeds<'a> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:891-900
    /// The seeds of `index` for the batch whose concatenated query is `query` (NCBI's
    /// `query->sequence`), with its contexts, their validity and their masks (the lookup
    /// segments: the unmasked intervals of each valid context in order, as NCBI's
    /// `BLAST_ComplementMaskLocations` gives them to the lookup table and to the index), and the
    /// word size (`lut_options->word_size`) and two-hit window (`word_options->window_size`, 0
    /// for the one-hit search) of `options`, and with `-db_hard_mask` the database and the
    /// algorithm of its masks (`hard_masks`).
    pub(super) fn new(
        index: &'a LoadedIndex,
        oids: OidMap<'a>,
        query: &'a [u8],
        contexts: &[QueryContext],
        context_valid: &[bool],
        options: SearchOptions,
        hard_masks: Option<(&'a super::database::BlastnDatabase, i32)>,
    ) -> Self {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
        // NCBI runs neither the preliminary search nor `DoPreSearch` for a batch (or a query
        // chunk) without a valid context.
        let search = if context_valid.iter().any(|&valid| valid) {
            index
                .searches
                .fetch_add(1, Ordering::Relaxed)
                .wrapping_add(1)
        } else {
            0
        };
        let visited: Vec<bool> = index
            .volumes
            .iter()
            .map(|volume| {
                oids.any_in(
                    volume.start_oid,
                    volume.start_oid.saturating_add(volume.n_oids),
                )
            })
            .collect();
        let plan = if search > 0 && !index.old_style {
            let mut held = index
                .held
                .lock()
                .expect("the held results are not poisoned");
            plan_search(&index.volumes, &visited, &mut held)
        } else {
            index
                .volumes
                .iter()
                .map(|_| VolumePlan::Searched { kept: false })
                .collect()
        };
        IndexSeeds {
            index,
            oids,
            query,
            segments: lookup_segments(contexts, context_valid),
            options,
            search,
            plan,
            visited,
            volumes: index.volumes.iter().map(|_| OnceLock::new()).collect(),
            moved: index
                .volumes
                .iter()
                .map(|_| AtomicBool::new(false))
                .collect(),
            scanned_hits: index.volumes.iter().map(|_| AtomicU32::new(0)).collect(),
            error: Mutex::new(None),
            hard_masks,
        }
    }

    /// Keeps the refusal `what` of volume `at` and kind `kind` unless one of an earlier volume,
    /// or of the same volume and a lower kind, is kept. Each kind's text depends only on the
    /// volume, so the refusal does not depend on the order of the search threads.
    fn fail(&self, at: usize, kind: u8, what: String) {
        let mut first = self.error.lock().expect("the error cell is not poisoned");
        if first.as_ref().is_none_or(|(key, _)| (at, kind) < *key) {
            *first = Some(((at, kind), what));
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:812-857
    /// The seeds of index volume `at`: those of an earlier search that it keeps, or the volume
    /// read and searched by the first caller (the others wait, as on NCBI's mutex) and unmapped
    /// after its search; `None` after an error, which is kept for [`Self::check`].
    fn volume(&self, at: usize) -> Option<&VolumeSeeds> {
        if let VolumePlan::Earlier(seeds) = &self.plan[at] {
            return Some(seeds);
        }
        self.volumes[at]
            .get_or_init(|| {
                let path = &self.index.volumes[at].path;
                let searched = IndexVolume::open(std::path::Path::new(path)).and_then(|volume| {
                    let masks = self.volume_masks(at);
                    search(&volume, self.query, &self.segments, self.options, &masks)
                        .map_err(|_| "a list or an offset outside its sections".to_string())
                });
                match searched {
                    Ok(seeds) => Some(Arc::new(seeds)),
                    Err(why) => {
                        self.fail(
                            at,
                            0,
                            format!(
                                "a MegaBLAST index volume that LOSAT does not read ({path}: {why})"
                            ),
                        );
                        None
                    }
                }
            })
            .as_deref()
    }

    /// With `-db_hard_mask`: the masked letters of the database sequences of index volume `at`
    /// that the search reads, by their sequence in the volume (`dbindex::search`), read with the
    /// masks of their blocks; a sequence whose masks do not read has none here (its own search
    /// reads them again and reports the error). The sequences that the search does not read (an
    /// OID list) have none. Out of line, as `dbindex::search`'s setup of the masks.
    #[inline(never)]
    fn volume_masks(&self, at: usize) -> Vec<(u32, Vec<(u32, u32)>)> {
        let Some((database, algo_id)) = self.hard_masks else {
            return Vec::new();
        };
        let volume = &self.index.volumes[at];
        let ordinals = self.oids.ordinals_in(
            volume.start_oid,
            volume.start_oid.saturating_add(volume.n_oids),
        );
        let mut masks = Vec::new();
        database.masks_of(ordinals, algo_id, |ordinal, stored| {
            if !stored.is_empty() {
                let seq = self.oids.oid(ordinal).get().wrapping_sub(volume.start_oid);
                masks.push((seq, stored.to_vec()));
            }
        });
        masks
    }

    /// The volume of database ordinal id `ordinal` and its sequence in the volume.
    fn locate(&self, ordinal: u32) -> (usize, u32) {
        self.index.locate(self.oids.oid(ordinal).get())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1727-1734
    /// Whether the database sequence `ordinal` is in an index volume (`check_oid` is not
    /// `eNotIndexed`): its chunks have the indexed word finder, else the scan.
    pub(super) fn is_indexed(&self, ordinal: u32) -> bool {
        self.index.volumes[self.locate(ordinal).0].has_index
    }

    /// Counts the `hits` initial hits that the scan of a chunk of the database sequence
    /// `ordinal` (of a volume without an index) saved, up to 2 for its volume.
    pub(super) fn note_scanned_hits(&self, ordinal: u32, hits: usize) {
        if hits == 0 {
            return;
        }
        let add = hits.min(2) as u32;
        let _ = self.scanned_hits[self.locate(ordinal).0].fetch_update(
            Ordering::Relaxed,
            Ordering::Relaxed,
            |count| Some(count.saturating_add(add).min(2)),
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:859-875
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1420-1423
    /// `CheckOid`: whether a chunk of the database sequence `ordinal` has seeds; NCBI neither
    /// reads nor searches a sequence without. A sequence of a volume without an index is searched
    /// (`eNotIndexed`).
    pub(super) fn has_seeds(&self, ordinal: u32) -> bool {
        let oid = self.oids.oid(ordinal).get();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:224-232
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:260-267
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:1059-1065
        // Each search of an old-style index appends the ends of its volumes to those of the
        // searches before and leaves their results null: in a later search an OID past the
        // index's sequences, below their number times the searches, finds a null result, and
        // NCBI's search fails ("Attempt to access NULL pointer").
        if self.index.old_style && self.search > 1 {
            let sequences: u64 = self
                .index
                .volumes
                .iter()
                .map(|volume| u64::from(volume.n_oids))
                .sum();
            if u64::from(oid) >= sequences
                && u64::from(oid) < sequences.saturating_mul(u64::from(self.search))
            {
                let first = &self.index.volumes[0].path;
                self.fail(
                    usize::MAX,
                    0,
                    format!(
                        "an old-style MegaBLAST index of fewer sequences than the database \
                         ({first}) and more than one query batch or query chunk"
                    ),
                );
                return false;
            }
        }
        let (at, subj) = self.index.locate(oid);
        if !self.index.volumes[at].has_index {
            return true;
        }
        self.volume(at).is_some_and(|seeds| seeds.has_seeds(subj))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:869-875
    /// `GetResults` of `oid` and `chunk`: the seeds of chunk `chunk` of the database sequence
    /// `ordinal`, the word size and the volume, or `None` (NCBI's word size 0); the word finder
    /// moves the list into its list of initial hits (noted for [`Self::check`]). A chunk that is
    /// a later sequence's in the index is refused (`VolumeSeeds::chunk_of_a_later_sequence`).
    fn chunk_seeds(&self, ordinal: u32, chunk: u32) -> Option<(&[Seed], u32, usize)> {
        let (at, subj) = self.locate(ordinal);
        let seeds = self.volume(at)?;
        if seeds.chunk_of_a_later_sequence(subj, chunk) {
            let path = &self.index.volumes[at].path;
            self.fail(
                at,
                1,
                format!(
                    "a MegaBLAST index volume that does not match the database ({path}: a \
                     sequence with more chunks than the index has)"
                ),
            );
            return None;
        }
        let list = seeds.seeds(subj, chunk)?;
        self.moved[at].store(true, Ordering::Relaxed);
        Some((list, seeds.word_size() as u32, at))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-262 (cited at
    // `core::subject_split::HardSplitState::next_chunk`)
    /// The seeds of the database sequence `ordinal` (of `length` letters) with the hard ranges
    /// `ranges` (`-db_hard_mask`, sorted; their chunks overlap by `chunk_overlap` within a range),
    /// placed into its hard chunks ([`HardSeeds`]; the Owner's answer C to Q-HARDIDX, approved
    /// exception 24 of PD-LOSAT-NCBI-DEFECTS): the lists of the index's chunks of the sequence's
    /// whole range (`GetResults` of each chunk number, as the search without the masks reads
    /// them, moved out as the word finder moves them), each seed at its letter of the sequence
    /// (its chunk's offset and its offset in the chunk). The refusals of the search without the
    /// masks are read here, for every chunk and seed: a chunk that is a later sequence's in the
    /// index, a seed at or past its chunk's last letter.
    pub(super) fn hard_seeds(
        &self,
        ordinal: u32,
        length: usize,
        ranges: &[(i32, i32)],
        chunk_overlap: usize,
    ) -> HardSeeds {
        let (at, subj) = self.locate(ordinal);
        let mut placed = HardSeeds::chunks_of(length, ranges, chunk_overlap, at);
        let Some(seeds) = self.volume(at) else {
            return placed;
        };
        placed.word_size = seeds.word_size() as u32;
        let mut whole = crate::core::subject_split::HardSplitState::default();
        whole.start(length, vec![(0, length as i32)], true);
        let mut number = 0u32;
        while let Some(chunk) = whole.next_chunk(chunk_overlap) {
            if seeds.chunk_of_a_later_sequence(subj, number) {
                let path = &self.index.volumes[at].path;
                self.fail(
                    at,
                    1,
                    format!(
                        "a MegaBLAST index volume that does not match the database ({path}: a \
                         sequence with more chunks than the index has)"
                    ),
                );
                return placed.cleared();
            }
            if let Some(list) = seeds.seeds(subj, number) {
                self.moved[at].store(true, Ordering::Relaxed);
                for seed in list {
                    if seed.s_off as usize >= chunk.letters() {
                        let path = &self.index.volumes[at].path;
                        self.fail(
                            at,
                            2,
                            format!(
                                "a MegaBLAST index seed past the last letter of its sequence \
                                 ({path})"
                            ),
                        );
                        return placed.cleared();
                    }
                    placed.place(*seed, chunk.offset + seed.s_off as usize, ranges);
                }
            }
            number += 1;
        }
        placed
    }

    /// After the preliminary stage: LOSAT's refusal of the search (the first by volume: an index
    /// volume that did not read or search, an index not made from the database, or a later
    /// search of an old-style index past its sequences; then a search whose scan saves initial
    /// hits into NCBI's list of no memory, [`Self::freed_list`]); else the results that the
    /// searched volumes keep for later searches, as the word finder left them.
    pub(super) fn check(&self) -> anyhow::Result<()> {
        // A volume that keeps its results is searched here if no sequence of it was (it is
        // reached), before the errors are read.
        for (at, plan) in self.plan.iter().enumerate() {
            if let VolumePlan::Searched { kept: true } = plan {
                let _ = self.volume(at);
            }
        }
        if let Some((_, what)) = self
            .error
            .lock()
            .expect("the error cell is not poisoned")
            .take()
        {
            return Err(super::mb_index::refusal(&what));
        }
        if self.freed_list() {
            let first = &self.index.volumes[0].path;
            return Err(super::mb_index::refusal(&format!(
                "a MegaBLAST database index of only some database volumes ({first}) and a later \
                 query batch or query chunk whose scan saves initial hits into a list that NCBI \
                 frees"
            )));
        }
        for at in 0..self.plan.len() {
            if let VolumePlan::Searched { kept: true } = self.plan[at] {
                let consumed = self.volume(at).map(|seeds| Arc::new(seeds.consumed()));
                let mut held = self
                    .index
                    .held
                    .lock()
                    .expect("the held results are not poisoned");
                held[at].results = consumed;
            }
        }
        Ok(())
    }

    /// Whether NCBI's search thread saves an initial hit into a freed list
    /// ([`saves_into_a_freed_list`] over the reached volumes).
    fn freed_list(&self) -> bool {
        saves_into_a_freed_list(
            self.plan
                .iter()
                .enumerate()
                .filter(|&(at, _)| self.visited[at])
                .map(|(at, plan)| {
                    (
                        plan,
                        self.moved[at].load(Ordering::Relaxed),
                        self.scanned_hits[at].load(Ordering::Relaxed),
                    )
                }),
        )
    }
}

/// One hard chunk of a database sequence ([`HardSeeds`]): its offset, residual and length (the
/// chunk loop's `SubjectChunk`, scanned from the residual), its hard range, and its seeds.
#[derive(Debug)]
struct HardSlot {
    offset: usize,
    residual: usize,
    length: usize,
    range: usize,
    seeds: Vec<Seed>,
}

/// The index's seeds of a database sequence of a search with `-db_hard_mask`, placed into its hard
/// chunks (the Owner's answer C of 2026-11-22 to Q-HARDIDX, 「Q5=C」: the seeds in their right
/// ranges, a result that NCBI BLAST+ 2.17.0 does not give, the same at every thread count;
/// approved exception 24 of PD-LOSAT-NCBI-DEFECTS). NCBI's word finder takes the list of the
/// index's chunk whose number is the hard chunk's (each hard range starts chunks of its own), so
/// its seeds land in other ranges or are taken from the next sequence's lists, by the order of its
/// threads. The search of the index found each seed in one unmasked stretch (`dbindex::search`:
/// a masked letter matches no query letter); here each seed goes to the hard range that holds its
/// left end, to the last chunk of the range that starts at or before it (the chunks of a range
/// over `MAX_DBSEQ_LEN` letters overlap; the index places a seed in the overlap in the later chunk
/// too), in the order of the lists and of each list. A sequence without masked letters thus has
/// NCBI's lists as the search without `-db_hard_mask` takes them.
#[derive(Debug)]
pub(super) struct HardSeeds {
    word_size: u32,
    at: usize,
    /// The chunks in the chunker's order (by range, then offset).
    chunks: Vec<HardSlot>,
}

impl HardSeeds {
    /// The hard chunks of a sequence of `length` letters with the hard ranges `ranges` (no seeds
    /// yet), of index volume `at`.
    fn chunks_of(
        length: usize,
        ranges: &[(i32, i32)],
        chunk_overlap: usize,
        at: usize,
    ) -> HardSeeds {
        let mut split = crate::core::subject_split::HardSplitState::default();
        split.start(length, ranges.to_vec(), true);
        let mut chunks = Vec::with_capacity(ranges.len());
        while let Some(chunk) = split.next_chunk(chunk_overlap) {
            chunks.push(HardSlot {
                offset: chunk.offset,
                residual: chunk.residual,
                length: chunk.letters(),
                range: chunk.range,
                seeds: Vec::new(),
            });
        }
        HardSeeds {
            word_size: 0,
            at,
            chunks,
        }
    }

    /// The chunks without seeds (after a refusal).
    fn cleared(mut self) -> HardSeeds {
        for chunk in &mut self.chunks {
            chunk.seeds.clear();
        }
        self
    }

    /// Places `seed`, whose left end is letter `start` of the sequence with the hard ranges
    /// `ranges` (sorted), into the chunks of the range that holds it. A seed in no range is not
    /// placed: the search finds none at a masked letter, so only a sequence whose masks did not
    /// read has one (its own search then fails).
    fn place(&mut self, seed: Seed, start: usize, ranges: &[(i32, i32)]) {
        let range = ranges.partition_point(|&(_, right)| right as usize <= start);
        if ranges
            .get(range)
            .is_some_and(|&(left, _)| left as usize <= start)
        {
            self.push(range, start, seed);
        }
    }

    /// Adds `seed`, at letter `start` of the sequence in hard range `range`, to the last chunk of
    /// the range that starts at or before it, with its offset in that chunk.
    fn push(&mut self, range: usize, start: usize, mut seed: Seed) {
        let after = self
            .chunks
            .partition_point(|chunk| (chunk.range, chunk.offset) <= (range, start));
        let chunk = &mut self.chunks[after - 1];
        debug_assert_eq!(chunk.range, range);
        seed.s_off = (start - chunk.offset) as u32;
        chunk.seeds.push(seed);
    }

    /// The seeds of the chunk at `offset` scanned from `residual`, of `length` letters, the word
    /// size and the volume; `None` for a chunk of none. Chunks of an equal offset, residual and
    /// length are those of hard ranges of no letters at one place, which hold no seed.
    fn chunk(
        &self,
        offset: usize,
        residual: usize,
        length: usize,
    ) -> Option<(&[Seed], u32, usize)> {
        let from = self.chunks.partition_point(|chunk| chunk.offset < offset);
        self.chunks[from..]
            .iter()
            .take_while(|chunk| chunk.offset == offset)
            .find(|chunk| (chunk.residual, chunk.length) == (residual, length))
            .filter(|chunk| !chunk.seeds.is_empty())
            .map(|chunk| (chunk.seeds.as_slice(), self.word_size, self.at))
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:248-259
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:325-356
/// Whether NCBI's search thread saves an initial hit into a freed list, given each reached
/// volume in order: what the search does with it, whether the word finder took a list of it, and
/// the initial hits (up to 2) that its scan saved. The thread's list of initial hits as the
/// volumes leave it: a list moved in from a volume searched here has memory; one from a volume
/// that keeps an earlier search's results has none (`allocated` 0); the first initial hit that a
/// scan then saves reallocates it to 0 bytes, which glibc gives one hit's block, and the next
/// one frees it.
fn saves_into_a_freed_list<'p>(volumes: impl Iterator<Item = (&'p VolumePlan, bool, u32)>) -> bool {
    #[derive(Clone, Copy)]
    enum List {
        Allocated,
        Empty,
        OneHit,
    }
    let mut list = List::Allocated;
    for (plan, moved, scanned_hits) in volumes {
        list = match (plan, list) {
            (VolumePlan::Searched { .. }, _) if moved => List::Allocated,
            (VolumePlan::Earlier(_), _) if moved => List::Empty,
            (VolumePlan::Scanned, List::Empty | List::OneHit) => match (list, scanned_hits) {
                (_, 0) => list,
                (List::Empty, 1) => List::OneHit,
                _ => return true,
            },
            _ => list,
        };
    }
    false
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1017-1058
/// The lookup segments (`[start, stop)` in the concatenated query): the unmasked intervals of
/// each valid context, in order (`BLAST_ComplementMaskLocations`; the contexts' masks are those
/// of the lookup tables, `lookup::build_unmasked_ranges`).
fn lookup_segments(contexts: &[QueryContext], context_valid: &[bool]) -> Vec<(u32, u32)> {
    let mut segments = Vec::new();
    for (ctx, &valid) in contexts.iter().zip(context_valid) {
        if !valid {
            continue;
        }
        let base = ctx.query_offset as u32;
        for (start, end) in build_unmasked_ranges(ctx.len, &ctx.masks) {
            segments.push((base + start as u32, base + end as u32));
        }
    }
    segments
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-262
/// NCBI's number of the subject chunk at `offset` (`subject->chunk`, counted from 0 by the chunk
/// loop, a chunk without ranges included): the chunks' offsets step by the chunk length less the
/// overlap, each rounded down to a whole byte of the packed letters. (With `-db_hard_mask` each
/// hard range starts chunks of its own, [`HardSeeds`].)
fn chunk_number(offset: usize, overlap: usize) -> u32 {
    let mut next = 0usize;
    let mut number = 0u32;
    loop {
        let start = next - next % COMPRESSION_RATIO;
        if start >= offset {
            return number;
        }
        next = start + MAX_DBSEQ_LEN - overlap;
        number += 1;
    }
}

/// The number of keys of the diagonal hash (`IR_HASH_SIZE`).
const IR_HASH_SIZE: u32 = 4 * 1024;

/// A diagonal of the hash (`ir_hash_entry`): its id, the query end of its last extension (0 for
/// none), and the next entry of its key (1-based in the pool, 0 for none).
#[derive(Clone, Copy, Default)]
struct IrEntry {
    diag: u32,
    qend: u32,
    next: u32,
}

/// The diagonal hash of one chunk (`ir_diag_hash`): an entry for each key, the entries added
/// after them, and the keys used (to clear them for the next chunk; NCBI makes a new hash for
/// each chunk).
struct IrHash {
    heads: Vec<IrEntry>,
    pool: Vec<IrEntry>,
    used: Vec<u32>,
}

impl IrHash {
    fn new() -> Self {
        IrHash {
            heads: vec![IrEntry::default(); IR_HASH_SIZE as usize],
            pool: Vec::new(),
            used: Vec::new(),
        }
    }

    /// The hash as NCBI makes it for a chunk: every entry empty.
    fn clear(&mut self) {
        for &key in &self.used {
            self.heads[key as usize] = IrEntry::default();
        }
        self.used.clear();
        self.pool.clear();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/index_ungapped.h:42-73
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/index_ungapped.c:120-153
    /// `IR_LOCATE` of `hash`, `diag` and `key`: the key's first entry when it is empty or holds
    /// `diag`; else the entry after it that holds `diag`, whose values it swaps with the first's,
    /// which it returns; else a new entry for `diag` after the first.
    fn locate(&mut self, diag: u32) -> &mut IrEntry {
        let key = diag % IR_HASH_SIZE;
        let head = self.heads[key as usize];
        if head.qend == 0 && head.diag == 0 && head.next == 0 {
            self.used.push(key);
        }
        if head.qend == 0 || head.diag == diag {
            return &mut self.heads[key as usize];
        }
        let mut at = head.next;
        while at != 0 {
            let entry = self.pool[at as usize - 1];
            if entry.diag == diag {
                let first = &mut self.heads[key as usize];
                let (diag, qend) = (first.diag, first.qend);
                first.diag = entry.diag;
                first.qend = entry.qend;
                let found = &mut self.pool[at as usize - 1];
                found.diag = diag;
                found.qend = qend;
                return &mut self.heads[key as usize];
            }
            at = entry.next;
        }
        self.pool.push(IrEntry {
            diag,
            qend: 0,
            next: head.next,
        });
        let at = self.pool.len() as u32;
        self.heads[key as usize].next = at;
        &mut self.pool[at as usize - 1]
    }
}

thread_local! {
    /// The diagonal hash of this thread's chunks, cleared after each.
    static IR_HASH: std::cell::RefCell<Option<IrHash>> = const { std::cell::RefCell::new(None) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1697-1782
/// `MB_IndexedWordFinder` of a chunk whose list holds `seeds` (in their order) for word size
/// `word_size`: the initial hits, with their ungapped extensions, into `hits`. The seeds of the
/// index run on past the end of a sequence into the next one of the index's local id, so a
/// seed's word can end past the chunk's last letter. `Err` for a seed that starts at or past the
/// chunk's last letter (in the last byte of its sequence, which the index's boundary codes keep
/// seeds from reaching in the searches LOSAT has seen, or past it, from an index not made from
/// the database): NCBI's later stages read such an extension past the sequence.
fn indexed_word_finder(
    scan: &ScanBatch<'_>,
    seeds: &[Seed],
    word_size: u32,
    s_seq_packed: &[u8],
    s_len: usize,
    hits: &mut Vec<UngappedHit>,
) -> Result<(), ()> {
    IR_HASH.with(|cell| {
        let mut cell = cell.borrow_mut();
        let hash = cell.get_or_insert_with(IrHash::new);
        let mut found = Ok(());
        for seed in seeds {
            let (q_off, s_off) = (seed.q_off, seed.s_off);
            // `IR_DIAG` of `q_off` and `s_off`.
            let diag = 0x1000_0000u32.wrapping_add(s_off).wrapping_sub(q_off);
            let entry = hash.locate(diag);
            if q_off.wrapping_add(word_size).wrapping_sub(1) <= entry.qend {
                continue;
            }
            if s_off as usize >= s_len {
                found = Err(());
                break;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-238
            let context = scan.query_context_index.context_for_offset(q_off as usize);
            let ctx = &scan.query_contexts[context];
            let ungapped = extend_hit_ungapped_approx_ncbi(
                scan.encoded_query_concat_blastna,
                s_seq_packed,
                q_off as usize,
                s_off as usize,
                s_off as usize + word_size as usize,
                s_len,
                scan.x_dropoff_scores[context],
                scan.nucl_score_table,
                scan.reduced_cutoff_scores[context],
                scan.score_matrix,
            );
            if ungapped.score >= scan.cutoff_scores[context] {
                x_push_ungapped(
                    hits,
                    UngappedHit::saved(
                        context as u32,
                        ctx,
                        q_off as usize,
                        s_off as usize,
                        &ungapped,
                    ),
                );
            }
            entry.diag = diag;
            // The `Int4` end of the extension, kept as `Uint4`.
            entry.qend = (ungapped.q_start as u32)
                .wrapping_add(ungapped.length as u32)
                .wrapping_sub(1);
        }
        hash.clear();
        found
    })
}

/// The search of one chunk of a database sequence with the index (`prelim_subject::
/// prelim_chunk_hits` with a used index): the chunk's seeds (none: no initial hit) through the
/// indexed word finder, then the chunk's HSP list as for every chunk, the saved initial hits not
/// counted. The diagonal table of the other word finders is not read, nor moved on.
#[inline(never)]
pub(super) fn index_chunk_hits<const UNGAPPED: bool>(
    scan: &ScanBatch<'_>,
    index: &IndexSeeds<'_>,
    subject: &PrelimSubject<'_>,
    chunk: &SubjectChunk,
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
    reuse_prelim_hits: bool,
) -> Vec<PrelimHit> {
    subject_scratch.prelim_hits.clear();
    subject_scratch.ungapped_hits.clear();
    let s_seq_packed = &subject.s_seq_packed_full[chunk.offset / COMPRESSION_RATIO..];
    let found = match subject.hard_index {
        // A hard chunk is scanned from its residual (`next_hard_chunk`).
        Some(placed) => {
            let residual = match chunk.seq_ranges {
                SeqRanges::Inline { range } => range.0 as usize,
                _ => 0,
            };
            placed.chunk(chunk.offset, residual, chunk.length)
        }
        None => index.chunk_seeds(
            subject.s_idx as u32,
            chunk_number(chunk.offset, subject.dbseq_chunk_overlap),
        ),
    };
    if let Some((seeds, word_size, at)) = found {
        if word_size > 0
            && indexed_word_finder(
                scan,
                seeds,
                word_size,
                s_seq_packed,
                chunk.length,
                &mut subject_scratch.ungapped_hits,
            )
            .is_err()
        {
            let path = &index.index.volumes[at].path;
            index.fail(
                at,
                2,
                format!("a MegaBLAST index seed past the last letter of its sequence ({path})"),
            );
            subject_scratch.ungapped_hits.clear();
        }
    }
    // NCBI's indexed word finder does not count the hits it saves (no
    // `Blast_UngappedStatsUpdate`).
    let uncounted = AtomicU64::new(0);
    chunk_hsp_list::<UNGAPPED>(
        scan.batch,
        subject,
        chunk,
        s_seq_packed,
        gap_scratch,
        subject_scratch,
        reuse_prelim_hits,
        &uncounted,
    )
}

#[cfg(test)]
#[path = "mb_index_seeds_tests.rs"]
mod tests;
