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
//! volume do not size the next query batch.

use std::sync::atomic::AtomicU64;
use std::sync::{Mutex, OnceLock};

use super::prelim_subject::{chunk_hsp_list, PrelimSubject, ScanBatch};
use super::{
    extend_hit_ungapped_approx_ncbi, x_push_ungapped, GapAlignScratch, PrelimHit, QueryContext,
    SubjectChunk, SubjectScratch, UngappedHit, COMPRESSION_RATIO, MAX_DBSEQ_LEN,
};
use crate::algorithm::blastn::lookup::build_unmasked_ranges;
use crate::dbindex::search::{search, LoadedIndex, SearchOptions, Seed, VolumeSeeds};
use crate::dbindex::volume::IndexVolume;
use crate::seqdb::OidMap;

/// The index's seeds of one query batch (`CIndexedDb_New` after `DoPreSearch`, or
/// `CIndexedDb_Old`): the index, the OIDs of the database's ordinal ids, the batch's
/// concatenated query and lookup segments, the word size, each volume's seeds once searched (or
/// `None` after an error), and the first error by volume.
pub(super) struct IndexSeeds<'a> {
    index: &'a LoadedIndex,
    oids: OidMap<'a>,
    query: &'a [u8],
    segments: Vec<(u32, u32)>,
    options: SearchOptions,
    volumes: Vec<OnceLock<Option<VolumeSeeds>>>,
    error: Mutex<Option<(usize, String)>>,
}

impl<'a> IndexSeeds<'a> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:891-900
    /// The seeds of `index` for the batch whose concatenated query is `query` (NCBI's
    /// `query->sequence`), with its contexts, their validity and their masks (the lookup
    /// segments: the unmasked intervals of each valid context in order, as NCBI's
    /// `BLAST_ComplementMaskLocations` gives them to the lookup table and to the index), and the
    /// word size (`lut_options->word_size`; the one-hit search: LOSAT refuses a two-hit one).
    pub(super) fn new(
        index: &'a LoadedIndex,
        oids: OidMap<'a>,
        query: &'a [u8],
        contexts: &[QueryContext],
        context_valid: &[bool],
        word_size: usize,
    ) -> Self {
        IndexSeeds {
            index,
            oids,
            query,
            segments: lookup_segments(contexts, context_valid),
            options: SearchOptions {
                word_size: word_size as u64,
                two_hits: 0,
            },
            volumes: index.volumes.iter().map(|_| OnceLock::new()).collect(),
            error: Mutex::new(None),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:812-857
    /// The seeds of volume `at`: the volume read and searched by the first caller (the others
    /// wait, as on NCBI's mutex) and unmapped after its search; `None` after an error, which is
    /// kept for [`Self::check`].
    fn volume(&self, at: usize) -> Option<&VolumeSeeds> {
        self.volumes[at]
            .get_or_init(|| {
                let path = &self.index.volumes[at].path;
                let searched = IndexVolume::open(std::path::Path::new(path)).and_then(|volume| {
                    search(&volume, self.query, &self.segments, self.options)
                        .map_err(|_| "a list or an offset outside its sections".to_string())
                });
                match searched {
                    Ok(seeds) => Some(seeds),
                    Err(why) => {
                        let mut first = self.error.lock().expect("the error cell is not poisoned");
                        if first.as_ref().is_none_or(|(volume, _)| at < *volume) {
                            *first = Some((at, format!("{path}: {why}")));
                        }
                        None
                    }
                }
            })
            .as_ref()
    }

    /// The volume of database ordinal id `ordinal` and its sequence in the volume.
    fn locate(&self, ordinal: u32) -> (usize, u32) {
        self.index.locate(self.oids.oid(ordinal).get())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:859-875
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1420-1423
    /// `CheckOid`: whether a chunk of the database sequence `ordinal` has seeds; NCBI neither
    /// reads nor searches a sequence without.
    pub(super) fn has_seeds(&self, ordinal: u32) -> bool {
        let (at, subj) = self.locate(ordinal);
        self.volume(at).is_some_and(|seeds| seeds.has_seeds(subj))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:869-875
    /// `GetResults(oid, chunk)`: the seeds of chunk `chunk` of the database sequence `ordinal` and
    /// the word size, or `None` (NCBI's word size 0).
    fn chunk_seeds(&self, ordinal: u32, chunk: u32) -> Option<(&[Seed], u32)> {
        let (at, subj) = self.locate(ordinal);
        let seeds = self.volume(at)?;
        Some((seeds.seeds(subj, chunk)?, seeds.word_size() as u32))
    }

    /// LOSAT's refusal of an index volume that did not read or search, after the preliminary
    /// stage (the first by volume).
    pub(super) fn check(&self) -> anyhow::Result<()> {
        match self
            .error
            .lock()
            .expect("the error cell is not poisoned")
            .take()
        {
            Some((_, what)) => Err(super::mb_index::refusal(&format!(
                "a MegaBLAST index volume that LOSAT does not read ({what})"
            ))),
            None => Ok(()),
        }
    }
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
/// overlap, each rounded down to a whole byte of the packed letters. (LOSAT refuses an index
/// with `-db_hard_mask`, whose ranges start chunks of their own.)
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
    /// `IR_LOCATE(hash, diag, key)`: the key's first entry when it is empty or holds `diag`;
    /// else the entry after it that holds `diag`, whose values it swaps with the first's, which
    /// it returns; else a new entry for `diag` after the first.
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
/// `word_size`: the initial hits, with their ungapped extensions, into `hits`.
fn indexed_word_finder(
    scan: &ScanBatch<'_>,
    seeds: &[Seed],
    word_size: u32,
    s_seq_packed: &[u8],
    s_len: usize,
    hits: &mut Vec<UngappedHit>,
) {
    IR_HASH.with(|cell| {
        let mut cell = cell.borrow_mut();
        let hash = cell.get_or_insert_with(IrHash::new);
        for seed in seeds {
            let (q_off, s_off) = (seed.q_off, seed.s_off);
            // `IR_DIAG(q_off, s_off)`.
            let diag = 0x1000_0000u32.wrapping_add(s_off).wrapping_sub(q_off);
            let entry = hash.locate(diag);
            if q_off.wrapping_add(word_size).wrapping_sub(1) <= entry.qend {
                continue;
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
    });
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
    let number = chunk_number(chunk.offset, subject.dbseq_chunk_overlap);
    if let Some((seeds, word_size)) = index.chunk_seeds(subject.s_idx as u32, number) {
        if word_size > 0 {
            indexed_word_finder(
                scan,
                seeds,
                word_size,
                s_seq_packed,
                chunk.length,
                &mut subject_scratch.ungapped_hits,
            );
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
