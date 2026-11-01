//! Backbone lookup table for TBLASTX - EXACT NCBI BLAST port
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c

use super::super::diagnostics::diagnostics_enabled;
use super::super::translation::QueryFrame;
use super::{
    compute_backbone_size, compute_mask, compute_unmasked_intervals, encode_kmer_3, ilog2, pv_set,
    QueryContext, LOOKUP_ALPHABET_SIZE, LOOKUP_WORD_LENGTH, PV_ARRAY_BTS,
};
use crate::config::ScoringMatrix;
use crate::stats::karlin_calc::{
    apply_check_ideal, compute_aa_composition, compute_karlin_params_ungapped,
    compute_score_freq_profile_for_matrix, compute_std_aa_composition,
};
use crate::stats::KarlinParams;
use crate::utils::matrix::{blosum62_ncbistdaa_score_row, protein_score};

pub const AA_HITS_PER_CELL: usize = 3;

// The BLASTX table filled in place, see x_lut_direct.rs (a child module, so that it uses this
// file's private builder pieces).
#[path = "x_lut_direct.rs"]
mod x_lut_direct;
// EXPERIMENT (LOSAT_X_TBXLUT): TBLASTX's full table built in place
// (x_lut_direct.rs); no NCBI counterpart for the re-export itself.
pub(crate) use x_lut_direct::x_tbxlut_reach_stub;
// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the spares of exact blastp's table built ahead of its block
// (x_lut_direct.rs); no NCBI counterpart for the re-export itself.
pub(crate) use x_lut_direct::XAaSpares;
// EXPERIMENT (LOSAT_X_BPLK16): the table of a blastp block, `u16` hits or `Int4` hits
// (x_lut_direct.rs); no NCBI counterpart for the re-export itself.
pub(crate) use x_lut_direct::XBpLookup;
// The comparison of two tables, for the unit tests (x_lut_direct.rs).
#[cfg(test)]
pub(crate) use x_lut_direct::{x_assert_same_table_as, x_same_table};
// The BLASTX table's vectors given back for the next build (LOSAT_X_BXPOOLN1:
// by searches of two threads or more).
pub(crate) use x_lut_direct::{x_bx_pooled_for, x_give_back_pooled};
// EXPERIMENT (PAIRLUT): the direct builder of the TBLASTX unit
// path's group tables (x_pair_lut.rs, a child module of x_lut_direct.rs), used by
// tblastx/x_pair_par.rs.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
// ```c
// BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
//                       &lookup_wrap->lut);
// ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
// BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
//                           query, lookup_segments, 0);
// ```
pub(crate) use x_lut_direct::x_pair_lut;

fn blosum62_ideal_karlin_params() -> KarlinParams {
    ideal_karlin_params_for_matrix(ScoringMatrix::Blosum62, (-4, 11))
}

// NCBI c++/src/algo/blast/core/blast_stat.c:1506-1528:
// for each finite matrix score, update sbp->loscore and sbp->hiscore.
fn matrix_score_bounds(matrix: ScoringMatrix) -> (i32, i32) {
    let mut low = i32::MAX;
    let mut high = i32::MIN;
    for query in 0..LOOKUP_ALPHABET_SIZE {
        for subject in 0..LOOKUP_ALPHABET_SIZE {
            let score = protein_score(matrix, query as u8, subject as u8);
            if score > i16::MIN as i32 && score < i16::MAX as i32 {
                low = low.min(score);
                high = high.max(score);
            }
        }
    }
    (low, high)
}

// NCBI c++/src/algo/blast/core/blast_stat.c:2833-2848:
// stdrfp = Blast_ResFreqNew(sbp); Blast_ResFreqStdComp(sbp, stdrfp);
// sfp = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
// BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
fn ideal_karlin_params_for_matrix(matrix: ScoringMatrix, bounds: (i32, i32)) -> KarlinParams {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2833-2848
    // ```c
    // stdrfp = Blast_ResFreqNew(sbp);
    // Blast_ResFreqStdComp(sbp, stdrfp);
    // sfp = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
    // BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
    // sbp->kbp_ideal = Blast_KarlinBlkNew();
    // Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
    // ```
    let std_comp = compute_std_aa_composition();
    let sfp =
        compute_score_freq_profile_for_matrix(&std_comp, &std_comp, bounds.0, bounds.1, matrix);
    compute_karlin_params_ungapped(&sfp)
        .expect("standard protein matrix ideal Karlin parameters must be computable")
}

/// The ideal Karlin block of `matrix` (`kbp_ideal`) with the score range that the profile
/// lookup's query preparation gives it (`x_lookup_table_for_profile`): BLOSUM62's -4..11, the
/// range of the matrix's scores otherwise. TBLASTN's composition redo takes its ungapped lambda
/// (`stage_d_kappa_params.rs`).
pub(crate) fn ideal_karlin_params_of(matrix: ScoringMatrix) -> KarlinParams {
    let bounds = if matrix == ScoringMatrix::Blosum62 {
        (-4, 11)
    } else {
        matrix_score_bounds(matrix)
    };
    ideal_karlin_params_for_matrix(matrix, bounds)
}

/// NCBI AaLookupBackboneCell - EXACT port
#[repr(C)]
#[derive(Clone, Copy)]
pub struct BackboneCell {
    pub num_used: i32,
    pub entries: [i32; AA_HITS_PER_CELL],
}

impl Default for BackboneCell {
    fn default() -> Self {
        Self {
            num_used: 0,
            entries: [0; AA_HITS_PER_CELL],
        }
    }
}

// No NCBI counterpart (LOSAT_X_BPLK16): NCBI's AaLookupBackboneCell holds Int4 offsets and counts.
// A blastp block's lookup query is at most about 50,000 letters (LOSAT_X_BPQBLOCK), so a query
// offset and a cell's hit count fit 16 bits and the cell takes 8 bytes instead of 16. The trait
// `XCell` is what the table's builder (x_lut_direct.rs) and the scan of the fused loop
// (blastp/x_seed_fused.rs) need of a cell; `BackboneCell` implements it with exactly the reads the
// code had before, so the tables of every other program stay what they were.
/// EXPERIMENT (LOSAT_X_BPLK16): a hit (a query offset) stored in a lookup cell or in the overflow
/// array: `i32` as NCBI's `Int4`, or `u16`.
pub trait XHit: Copy + Default + Send + Sync + 'static {
    /// The query offset.
    fn q(self) -> i32;

    /// `dst` set to the query offsets `src` (`dst.len() == src.len()`); a `u16` hit panics on an
    /// offset above 65,535.
    fn x_store(dst: &mut [Self], src: &[i32]);

    /// The hit at `p` read with one volatile load (not taken out of a vector register by the
    /// compiler, `x_cell_avx2_defer_masked`).
    ///
    /// # Safety
    /// `p` must point to a hit.
    unsafe fn x_read_volatile(p: *const Self) -> i32;

    /// Eight hits from `p` as eight `i32` lanes: one `vmovdqu` for `i32`, one `vpmovzxwd ymm, m128`
    /// for `u16`.
    ///
    /// # Safety
    /// The CPU must have AVX2 and `p` must point to eight hits.
    #[cfg(target_arch = "x86_64")]
    unsafe fn x_load8(p: *const Self) -> std::arch::x86_64::__m256i;

    /// The lanes of `x_load8` for the hits at `p` where `in_rest` is set (all ones per lane); the
    /// other lanes are not used by the caller. `i32`: `vpmaskmovd`, which reads nothing past the
    /// hits. `u16`: the same 16-byte widening load as `x_load8`, so the 16 bytes at `p` must lie
    /// inside the allocation (`XCell::PAD_HITS`, `XCell::PAD_CELLS`).
    ///
    /// # Safety
    /// The CPU must have AVX2, `in_rest` must be set exactly for the first `n < 8` lanes and `p`
    /// must point to `n` hits, with (for `u16`) 16 readable bytes from `p` on.
    #[cfg(target_arch = "x86_64")]
    unsafe fn x_load_tail(
        p: *const Self,
        in_rest: std::arch::x86_64::__m256i,
    ) -> std::arch::x86_64::__m256i;
}

impl XHit for i32 {
    #[inline(always)]
    fn q(self) -> i32 {
        self
    }

    #[inline(always)]
    fn x_store(dst: &mut [Self], src: &[i32]) {
        dst.copy_from_slice(src);
    }

    #[inline(always)]
    unsafe fn x_read_volatile(p: *const Self) -> i32 {
        // SAFETY: (caller) `p` points to a hit.
        unsafe { std::ptr::read_volatile(p) }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn x_load8(p: *const Self) -> std::arch::x86_64::__m256i {
        // SAFETY: (caller) the CPU has AVX2 and `p` points to eight hits.
        unsafe { std::arch::x86_64::_mm256_loadu_si256(p.cast()) }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn x_load_tail(
        p: *const Self,
        in_rest: std::arch::x86_64::__m256i,
    ) -> std::arch::x86_64::__m256i {
        // SAFETY: (caller) the CPU has AVX2; the masked load reads only the lanes in `in_rest`.
        unsafe { std::arch::x86_64::_mm256_maskload_epi32(p, in_rest) }
    }
}

impl XHit for u16 {
    #[inline(always)]
    fn q(self) -> i32 {
        i32::from(self)
    }

    fn x_store(dst: &mut [Self], src: &[i32]) {
        for (d, &s) in dst.iter_mut().zip(src) {
            *d = u16::try_from(s).expect("a query offset of a u16 lookup table fits 16 bits");
        }
    }

    #[inline(always)]
    unsafe fn x_read_volatile(p: *const Self) -> i32 {
        // SAFETY: (caller) `p` points to a hit.
        i32::from(unsafe { std::ptr::read_volatile(p) })
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn x_load8(p: *const Self) -> std::arch::x86_64::__m256i {
        use std::arch::x86_64::*;
        // SAFETY: (caller) the CPU has AVX2 and `p` points to eight hits (16 bytes).
        unsafe { _mm256_cvtepu16_epi32(_mm_loadu_si128(p.cast())) }
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn x_load_tail(
        p: *const Self,
        _in_rest: std::arch::x86_64::__m256i,
    ) -> std::arch::x86_64::__m256i {
        use std::arch::x86_64::*;
        // SAFETY: (caller) the CPU has AVX2 and the 16 bytes at `p` are inside the allocation;
        // the lanes past the hits are blended away by the caller.
        unsafe { _mm256_cvtepu16_epi32(_mm_loadu_si128(p.cast())) }
    }
}

/// EXPERIMENT (LOSAT_X_BPLK16): a cell of the lookup table's backbone: `BackboneCell` (`Int4`
/// hits, as NCBI's AaLookupBackboneCell) or `BackboneCell16`. A cell with more than
/// `AA_HITS_PER_CELL` hits keeps the cursor of its hits in the table's overflow array in its
/// entries.
pub trait XCell: Copy + Default + Send + Sync + 'static {
    /// The cell's hits.
    type Hit: XHit;

    /// Zero cells kept after the table's last cell (the AVX2 tail of the scan reads 16 bytes from
    /// a hit list's start, which for a cell's own entries runs into the cells after it).
    const PAD_CELLS: usize;

    /// Zero hits kept after the overflow array's last hit (the same read, at the end of the
    /// overflow array).
    const PAD_HITS: usize;

    /// The cell's `num_used`.
    fn count(&self) -> usize;

    /// The hits and the used cells (`count()` above 0) of `cells`.
    fn x_hits_and_used(cells: &[Self]) -> (u64, u64) {
        let (mut hits, mut used) = (0u64, 0u64);
        for cell in cells {
            if cell.count() > 0 {
                hits += cell.count() as u64;
                used += 1;
            }
        }
        (hits, used)
    }

    /// The cell's hits: its own entries, or the `count()` hits from its cursor in `overflow`.
    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:60-99 (quoted at `x_scan_cells`,
    // blastp/x_seed_fused.rs)
    fn hits<'a>(&'a self, overflow: &'a [Self::Hit]) -> &'a [Self::Hit];

    /// Sets `num_used`.
    fn set_count(&mut self, count: usize);

    /// Sets the cursor of a cell with more than `AA_HITS_PER_CELL` hits.
    fn set_cursor(&mut self, cursor: usize);

    /// The cell's own entries (a cell with at most `AA_HITS_PER_CELL` hits is written here).
    fn entries_mut(&mut self) -> &mut [Self::Hit];
}

impl XCell for BackboneCell {
    type Hit = i32;
    const PAD_CELLS: usize = 0;
    const PAD_HITS: usize = 0;

    #[inline(always)]
    fn count(&self) -> usize {
        self.num_used as usize
    }

    // The loop `x_bp_table_ratio` had before LOSAT_X_BPLK16, for the instance of the table of
    // 32-bit hits to keep its code.
    fn x_hits_and_used(cells: &[Self]) -> (u64, u64) {
        let (mut hits, mut used) = (0u64, 0u64);
        for cell in cells {
            if cell.num_used > 0 {
                hits += cell.num_used as u64;
                used += 1;
            }
        }
        (hits, used)
    }

    #[inline(always)]
    fn hits<'a>(&'a self, overflow: &'a [i32]) -> &'a [i32] {
        let numhits = self.num_used as usize;
        if numhits <= AA_HITS_PER_CELL {
            &self.entries[..numhits]
        } else {
            let cursor = self.entries[0] as usize;
            &overflow[cursor..cursor + numhits]
        }
    }

    #[inline(always)]
    fn set_count(&mut self, count: usize) {
        self.num_used = count as i32;
    }

    #[inline(always)]
    fn set_cursor(&mut self, cursor: usize) {
        self.entries[0] = cursor as i32;
    }

    #[inline(always)]
    fn entries_mut(&mut self) -> &mut [i32] {
        &mut self.entries
    }
}

/// EXPERIMENT (LOSAT_X_BPLK16): the cell of the lookup table of a blastp block whose lookup
/// query is at most 65,535 letters (`LOSAT_X_BPLK16`): `BackboneCell` with 16-bit hits and count,
/// 8 bytes. A cell with more than `AA_HITS_PER_CELL` hits keeps its cursor in the overflow array as
/// a `u32` in `entries[0]` (low half) and `entries[1]` (high half).
#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct BackboneCell16 {
    pub num_used: u16,
    pub entries: [u16; AA_HITS_PER_CELL],
}

const _: () = assert!(std::mem::size_of::<BackboneCell16>() == 8);

impl XCell for BackboneCell16 {
    type Hit = u16;
    // A 16-byte load at a cell's entries (2 bytes into the cell) ends 10 bytes past the cell's
    // end at most.
    const PAD_CELLS: usize = 2;
    // A 16-byte load at the start of an overflow hit list of at least one hit ends 7 hits past
    // its end at most.
    const PAD_HITS: usize = 8;

    #[inline(always)]
    fn count(&self) -> usize {
        usize::from(self.num_used)
    }

    #[inline(always)]
    fn hits<'a>(&'a self, overflow: &'a [u16]) -> &'a [u16] {
        let numhits = usize::from(self.num_used);
        if numhits <= AA_HITS_PER_CELL {
            &self.entries[..numhits]
        } else {
            let cursor = usize::from(self.entries[0]) | usize::from(self.entries[1]) << 16;
            &overflow[cursor..cursor + numhits]
        }
    }

    #[inline(always)]
    fn set_count(&mut self, count: usize) {
        self.num_used = u16::try_from(count).expect("a cell of a u16 lookup table counts 16 bits");
    }

    #[inline(always)]
    fn set_cursor(&mut self, cursor: usize) {
        let cursor =
            u32::try_from(cursor).expect("the overflow cursor of a u16 table fits 32 bits");
        self.entries[0] = cursor as u16;
        self.entries[1] = (cursor >> 16) as u16;
    }

    #[inline(always)]
    fn entries_mut(&mut self) -> &mut [u16] {
        &mut self.entries
    }
}

// NCBI c++/src/algo/blast/core/blast_aalookup.c:562-603:
// row = info->matrix[query_word[current_pos]];
// score + row[i] controls neighboring word inclusion.
#[inline(always)]
fn lookup_matrix_score(matrix: ScoringMatrix, query: usize, subject: usize) -> i32 {
    if matrix == ScoringMatrix::Blosum62 {
        blosum62_ncbistdaa_score_row(query)[subject]
    } else {
        protein_score(matrix, query as u8, subject as u8)
    }
}

/// NCBI-style lookup table with flat backbone
///
/// EXPERIMENT (LOSAT_X_BPLK16): generic over its cell, `BackboneCell` (the default: every table of
/// every program) or `BackboneCell16` (blastp's block tables with `u16` hits). A `BackboneCell16`
/// table has `XCell::PAD_CELLS` zero cells after its `backbone_size` cells and `XCell::PAD_HITS`
/// zero hits after its overflow (the AVX2 tail of the fused loop reads 16 bytes from a hit list's
/// start); the padding is in no count, `longest_chain` or presence bit.
pub struct BlastAaLookupTable<C: XCell = BackboneCell> {
    pub backbone: Vec<C>,
    pub overflow: Vec<C::Hit>,
    // Presence vector: PV_ARRAY_TYPE = Uint4
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:41-44
    pub pv: Vec<u32>,
    pub frame_bases: Vec<i32>,
    pub num_contexts: usize,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:246-250
    // ```c
    // Uint4 QueryInfo_GetSeqBufLen(const BlastQueryInfo* qinfo)
    // {
    //     BlastContextInfo * cinfo = & qinfo->contexts[qinfo->last_context];
    //     return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-127
    // ```c
    // (*seq_blk)->sequence_start = (Uint1 *) buffer;
    // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
    // (*seq_blk)->length = length;
    // ```
    pub query_length: i32,
    pub word_length: i32,
    pub alphabet_size: i32,
    pub charsize: i32,
    pub mask: i32,
    pub longest_chain: i32,
    // NCBI BlastAaLookupTable has no lazy-neighbor flag; only NCBI fields are stored.
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:95-140
    pub threshold: i32,
    pub row_max: Vec<i32>,
}

impl<C: XCell> BlastAaLookupTable<C> {
    #[inline]
    pub fn get_context_idx(&self, concat_off: i32) -> usize {
        let bases = &self.frame_bases;
        let mut lo = 0usize;
        let mut hi = self.num_contexts;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if concat_off < bases[mid] {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        lo.saturating_sub(1)
    }
}

impl BlastAaLookupTable {
    #[inline(always)]
    pub fn get_hits(&self, index: usize) -> &[i32] {
        let cell = unsafe { self.backbone.get_unchecked(index) };
        let num = cell.num_used as usize;
        if num == 0 {
            return &[];
        }
        if num <= AA_HITS_PER_CELL {
            unsafe { std::slice::from_raw_parts(cell.entries.as_ptr(), num) }
        } else {
            let cursor = cell.entries[0] as usize;
            unsafe { std::slice::from_raw_parts(self.overflow.as_ptr().add(cursor), num) }
        }
    }
}

// The per-cell chains of a thin backbone are carved from a few large blocks owned by the backbone
// instead of being one `malloc`/`realloc`/`free` each.
//
// A lookup table for one query chunk makes tens of thousands of chains that start at 8 `Int4` and
// double as they fill, and frees them all when the final backbone has been written. Under
// WebAssembly threads the allocator has one lock for all threads, and two threads building lookup
// tables spend more time waiting for it than allocating. The arena hands out blocks of the same
// sizes with the same contents (zeroed where `calloc`/`Vec::resize` zero them) and reuses the
// block a chain leaves when it doubles, so every chain holds exactly the integers a `Vec` chain
// held; only where the block lives changes. Everything is released when the backbone is dropped.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:48-72
// ```c
// if (backbone[index] == NULL) {
//     chain_size = 8;
//     hits_in_chain = 0;
//     chain = (Int4 *) malloc(chain_size * sizeof(Int4));
// ...
//     backbone[index] = chain;
// ...
// if ((hits_in_chain + 2) == chain_size) {
//     chain_size = chain_size * 2;
//     chain = (Int4 *) realloc(chain, chain_size * sizeof(Int4));
// ...
//     backbone[index] = chain;
//     chain[0] = chain_size;
// ```
// NCBI makes one malloc for each new chain and one realloc for each doubling. The arena changes
// where those blocks live and nothing else: a new chain still has 8 Int4, a full chain still
// doubles, and every chain holds the same integers in the same order (placement only).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:356-360
// ```c
// for (j=0; j <lookup->thin_backbone[i][1]; j++)
//     dest[j] = lookup->thin_backbone[i][j + 2];
// /* done with this chain- free it */
// sfree(lookup->thin_backbone[i]);
// lookup->thin_backbone[i] = NULL;
// ```
// NCBI frees each chain when it has been copied into the final backbone. The arena frees all its
// blocks together when the backbone is dropped. Reuse of memory only.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:52,67
// ```c
// chain = (Int4 *) malloc(chain_size * sizeof(Int4));
// ...
// chain = (Int4 *) realloc(chain, chain_size * sizeof(Int4));
// ```
// The blocks of the arena take the place of the malloc and realloc of the C code above. Their size
// classes are 8, 16, 32, ... Int4, the sizes the chain has in C.
/// Power-of-two blocks of `Int4` for the chains of one backbone.
struct XChainArena {
    /// Every allocation made for this arena; freed with it.
    slabs: Vec<Box<[i32]>>,
    cursor: *mut i32,
    remaining: usize,
    /// Blocks chains have left, by size class (`1 << class` ints).
    free: Vec<Vec<std::ptr::NonNull<i32>>>,
    /// Slabs of `SLAB_INTS` kept by `x_reset` for the next build (the pooled build).
    spare: Vec<Box<[i32]>>,
    /// EXPERIMENT (LOSAT_X_BPDBKEEP): every allocation is made by `x_allocate_unthresholded`.
    x_unthresholded: bool,
    /// EXPERIMENT (LOSAT_X_BPDBKEEP): the regions of `x_allocate_unthresholded`.
    x_unthresholded_slabs: Vec<Vec<i32>>,
    /// EXPERIMENT (LOSAT_X_BPDBKEEP): the ints of the last region given out.
    x_unthresholded_used: usize,
}

impl XChainArena {
    /// `Int4`s per slab (256 KiB).
    const SLAB_INTS: usize = 1 << 16;

    fn new() -> Self {
        Self {
            slabs: Vec::new(),
            cursor: std::ptr::null_mut(),
            remaining: 0,
            free: Vec::new(),
            spare: Vec::new(),
            x_unthresholded: false,
            x_unthresholded_slabs: Vec::new(),
            x_unthresholded_used: 0,
        }
    }

    // No NCBI counterpart (buffer reuse, the pooled build): every block given up at once, the
    // slabs kept.
    fn x_reset(&mut self) {
        for slab in self.slabs.drain(..) {
            if slab.len() == Self::SLAB_INTS {
                self.spare.push(slab);
            }
        }
        self.x_unthresholded_slabs.clear();
        self.x_unthresholded_used = 0;
        self.cursor = std::ptr::null_mut();
        self.remaining = 0;
        for list in &mut self.free {
            list.clear();
        }
    }

    fn allocate(&mut self, ints: usize) -> *mut i32 {
        if ints == Self::SLAB_INTS {
            if let Some(mut slab) = self.spare.pop() {
                let start = slab.as_mut_ptr();
                self.slabs.push(slab);
                return start;
            }
        }
        if self.x_unthresholded {
            return self.x_allocate_unthresholded(ints);
        }
        let mut slab = vec![0i32; ints].into_boxed_slice();
        let start = slab.as_mut_ptr();
        // Moving the box into the list does not move the ints.
        self.slabs.push(slab);
        start
    }

    // No NCBI counterpart (LOSAT_X_BPDBKEEP): where the ints are allocated.
    /// `allocate(ints)` from a region that glibc maps and unmaps without raising its mmap threshold
    /// (`utils::x_mapped::x_unthresholded_capacity`: zeroed by `calloc`, whose mapped pages are not
    /// touched until written). The slabs and long chains are cut from the last region while it has
    /// room, so that a build reserves the address space of one region per 32 MiB of chains, not
    /// one per slab (3,000 human proteins x `BlastDBSpecies6` at one thread, blocks of 50,000 query
    /// letters: VmPeak 129 -> 97 MB, RssAnon the same).
    #[inline(never)]
    fn x_allocate_unthresholded(&mut self, ints: usize) -> *mut i32 {
        let room = self
            .x_unthresholded_slabs
            .last()
            .map_or(0, |region| region.len() - self.x_unthresholded_used);
        if room < ints {
            self.x_unthresholded_slabs.push(vec![
                0i32;
                crate::utils::x_mapped::x_unthresholded_capacity(
                    ints,
                    std::mem::size_of::<i32>()
                )
            ]);
            self.x_unthresholded_used = 0;
        }
        let used = self.x_unthresholded_used;
        let region = self
            .x_unthresholded_slabs
            .last_mut()
            .expect("a region with room for the ints");
        debug_assert!(used + ints <= region.len());
        // SAFETY: `used + ints <= region.len()`, and the region is never resized, so the ints
        // stay inside it, where earlier cuts do not reach, until the arena drops it.
        let start = unsafe { region.as_mut_ptr().add(used) };
        self.x_unthresholded_used = used + ints;
        start
    }

    /// A block of `1 << class` ints with unspecified contents.
    fn take(&mut self, class: u8) -> std::ptr::NonNull<i32> {
        let ints = 1usize << class;
        if let Some(block) = self
            .free
            .get_mut(class as usize)
            .and_then(|list| list.pop())
        {
            return block;
        }
        let block = if ints > Self::SLAB_INTS / 4 {
            self.allocate(ints)
        } else {
            if self.remaining < ints {
                self.cursor = self.allocate(Self::SLAB_INTS);
                self.remaining = Self::SLAB_INTS;
            }
            let block = self.cursor;
            // SAFETY: `ints <= remaining`, so the block stays inside the slab.
            self.cursor = unsafe { block.add(ints) };
            self.remaining -= ints;
            block
        };
        std::ptr::NonNull::new(block).expect("allocation is not null")
    }

    /// Return a block obtained from `take(class)` of this arena.
    fn give(&mut self, block: std::ptr::NonNull<i32>, class: u8) {
        if self.free.len() <= class as usize {
            self.free.resize_with(class as usize + 1, Vec::new);
        }
        self.free[class as usize].push(block);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:49-56
// ```c
// if (backbone[index] == NULL) {
//     chain_size = 8;
//     hits_in_chain = 0;
//     chain = (Int4 *) malloc(chain_size * sizeof(Int4));
//     ASSERT(chain != NULL);
//     chain[0] = chain_size;
//     chain[1] = hits_in_chain;
//     backbone[index] = chain;
// ```
// A chain is the C array: chain[0] = allocated size, chain[1] = number of hits, then the hits.
// It is a block of the arena of the backbone that holds it.
/// One chain of a thin backbone: `chain[0]` = allocated ints, `chain[1]` =
/// hits, then the hits (NCBI `Int4 *`): `1 << class` ints inside the arena of
/// the backbone that holds this chain. Only ever stored in that backbone's
/// `cells`.
struct LookupChain {
    block: std::ptr::NonNull<i32>,
    class: u8,
}

impl LookupChain {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:49-56
    // ```c
    // chain_size = 8;
    // hits_in_chain = 0;
    // chain = (Int4 *) malloc(chain_size * sizeof(Int4));
    // ASSERT(chain != NULL);
    // chain[0] = chain_size;
    // chain[1] = hits_in_chain;
    // backbone[index] = chain;
    // ```
    // The new chain comes from the arena: chain_size 8, all ints zero.
    /// The 8-int chain of blast_lookup.c:50-56. NCBI `malloc`s it and sets `chain[0..2]`; the
    /// arena block is zeroed in all 8 ints.
    fn x_zeroed8(arena: &mut XChainArena) -> Self {
        let block = arena.take(3);
        // SAFETY: the block has 8 ints.
        unsafe { std::ptr::write_bytes(block.as_ptr(), 0, 8) };
        LookupChain { block, class: 3 }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:65-72
    // ```c
    // if ((hits_in_chain + 2) == chain_size) {
    //     chain_size = chain_size * 2;
    //     chain = (Int4 *) realloc(chain, chain_size * sizeof(Int4));
    // ...
    //     backbone[index] = chain;
    //     chain[0] = chain_size;
    // ```
    // A chain is copied into a block of twice the size and the old block is given back. The
    // first old_size ints are the same, the rest are zero.
    /// `Vec::resize(new_size, 0)` for a chain that doubles.
    fn x_resize(&mut self, new_size: usize, arena: &mut XChainArena) {
        let old_size = 1usize << self.class;
        assert!(
            new_size == old_size * 2,
            "lookup chains double: {old_size} -> {new_size}"
        );
        let new_block = arena.take(self.class + 1);
        // SAFETY: both blocks are valid for their sizes and distinct.
        unsafe {
            std::ptr::copy_nonoverlapping(self.block.as_ptr(), new_block.as_ptr(), old_size);
            std::ptr::write_bytes(new_block.as_ptr().add(old_size), 0, old_size);
        }
        arena.give(self.block, self.class);
        self.block = new_block;
        self.class += 1;
    }
}

impl std::ops::Deref for LookupChain {
    type Target = [i32];
    #[inline(always)]
    fn deref(&self) -> &[i32] {
        // SAFETY: the block has `1 << class` initialised ints, and the arena it lies in lives as
        // long as the backbone this chain is borrowed from.
        unsafe { std::slice::from_raw_parts(self.block.as_ptr(), 1usize << self.class) }
    }
}

impl std::ops::DerefMut for LookupChain {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut [i32] {
        // SAFETY: as in `deref`, and `&mut self` is the only access.
        unsafe { std::slice::from_raw_parts_mut(self.block.as_ptr(), 1usize << self.class) }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:246-247
// ```c
// lookup->thin_backbone =
//     (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
// ```
// The backbone of chains is the calloc of `backbone_size` pointers in the C code (here all None).
/// The chains of one thin backbone (NCBI `Int4 **`), with the memory of the
/// pooled ones. `cells` is declared first so that it is dropped first.
struct LookupBackboneChains {
    cells: Vec<Option<LookupChain>>,
    arena: XChainArena,
}

impl LookupBackboneChains {
    /// `calloc(backbone_size, sizeof(Int4 *))`
    fn new(backbone_size: usize) -> Self {
        Self {
            cells: std::iter::repeat_with(|| None)
                .take(backbone_size)
                .collect(),
            arena: XChainArena::new(),
        }
    }

    // No NCBI counterpart (LOSAT_X_BPDBKEEP): where the cells and the chains are allocated.
    /// `new(backbone_size)` whose cells and slabs glibc maps and unmaps without raising its mmap
    /// threshold (`utils::x_mapped::x_unthresholded_capacity`): a build's one-off 446 KiB of
    /// cells (BLASTP's 28,540) and 256 KiB slabs would raise it to their size when they are freed.
    fn x_new_unthresholded(backbone_size: usize) -> Self {
        let mut cells = Vec::with_capacity(crate::utils::x_mapped::x_unthresholded_capacity(
            backbone_size,
            std::mem::size_of::<Option<LookupChain>>(),
        ));
        cells.resize_with(backbone_size, || None);
        let mut arena = XChainArena::new();
        arena.x_unthresholded = true;
        Self { cells, arena }
    }

    // No NCBI counterpart (buffer reuse, the pooled build): as `new(backbone_size)`, on this
    // backbone's vectors and slabs.
    fn x_reset(&mut self, backbone_size: usize) {
        self.cells.clear();
        self.cells.resize_with(backbone_size, || None);
        self.arena.x_reset();
    }
}

impl std::ops::Deref for LookupBackboneChains {
    type Target = [Option<LookupChain>];
    #[inline(always)]
    fn deref(&self) -> &[Option<LookupChain>] {
        &self.cells
    }
}

#[derive(Default)]
struct LookupBuildStats {
    total_exact_positions: usize,
    skipped_invalid_residue: usize,
    skipped_seg_mask: usize,
    exact_added_count: usize,
    neighbor_added_count: usize,
    words_with_exact_only: usize,
    words_with_neighbors: usize,
    max_neighbors_for_single_word: usize,
    max_neighbor_word_idx: usize,
    neighbor_words_generated: usize,
}

struct NeighborInfo<'a> {
    thin_backbone: &'a mut LookupBackboneChains,
    query_word: &'a [u8],
    subject_word: [u8; 32],
    alphabet_size: usize,
    wordsize: usize,
    charsize: usize,
    row_max: &'a [i32],
    offset_list: &'a [i32],
    threshold: i32,
    query_bias: i32,
    matrix: ScoringMatrix,
    neighbor_added_count: usize,
    neighbor_words_generated: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:96-97
// ```c
// BlastAaLookupIndexQuery((BlastAaLookupTable*) lookup_wrap->lut, matrix,
//                         query, lookup_segments, 0);
// ```
struct PreparedLookupQuery {
    concat_query: Vec<u8>,
    lookup_locations: Vec<(i32, i32)>,
    frame_bases: Vec<i32>,
    contexts: Vec<QueryContext>,
    skipped_seg_mask: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1032-1114
// ```c
// start_offset = query_info->contexts[context].query_offset;
// end_offset = query_info->contexts[context].query_length + start_offset - 1;
// ...
// if (mask_loc == NULL || mask_loc->seqloc_array[context] == NULL) {
//     tail = BlastSeqLocNew(tail ? &tail : complement_mask,
//                           start_offset, end_offset);
//     continue;
// }
// ...
// filter_start = start_offset + seq_range->left;
// filter_end = start_offset + seq_range->right;
// ...
// tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
// ```
fn build_lookup_segments_for_context(
    query_offset: i32,
    query_length: usize,
    seg_masks: &[(usize, usize)],
) -> Vec<(i32, i32)> {
    if query_length == 0 {
        return Vec::new();
    }

    let start_offset = query_offset;
    let end_offset = query_offset
        + i32::try_from(query_length).expect("NCBI BLAST query length must fit in Int4")
        - 1;

    if seg_masks.is_empty() {
        return vec![(start_offset, end_offset)];
    }

    let mut sorted_masks = seg_masks.to_vec();
    sorted_masks.sort_by_key(|&(start, _)| start);

    let mut lookup_segments = Vec::new();
    let mut first = true;
    let mut last_interval_open = true;
    let mut left = 0i32;

    for (start, end) in sorted_masks {
        let filter_start =
            start_offset + i32::try_from(start).expect("NCBI BLAST mask start must fit in Int4");
        let filter_end =
            start_offset + i32::try_from(end).expect("NCBI BLAST mask end must fit in Int4") - 1;

        if first {
            last_interval_open = true;
            first = false;

            if filter_start > start_offset {
                left = start_offset;
            } else {
                left = filter_end + 1;
                continue;
            }
        }

        let right = filter_start - 1;
        if left <= right {
            lookup_segments.push((left, right));
        }
        if filter_end >= end_offset {
            last_interval_open = false;
            break;
        }
        left = filter_end + 1;
    }

    if last_interval_open && left <= end_offset {
        lookup_segments.push((left, end_offset));
    }

    lookup_segments
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:96-100
// ```c
// if (word_length > to - from + 1)
//     continue;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:102-128
// ```c
// seq = query->sequence + from;
// word_target = seq + lut_word_length;
// for (offset = from; offset <= to; offset++, seq++) { ... }
// ```
#[inline]
fn count_lookup_positions(lookup_segments: &[(i32, i32)], word_length: usize) -> usize {
    lookup_segments
        .iter()
        .map(|&(left, right)| {
            let seg_len = usize::try_from(right - left + 1)
                .expect("NCBI BLAST lookup segment length must be non-negative");
            seg_len.saturating_sub(word_length - 1)
        })
        .sum()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-314
// ```c
// info->contexts[i].query_offset = new_offsets[i];
// distance = new_offsets[i+1] - new_offsets[i];
// info->contexts[i].query_length = distance ? distance-1 : 0;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:621-636
// ```c
// if (!mask_at_hash) {
//     BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc, program_number);
// }
// if (lookup_segments) {
//     BLAST_ComplementMaskLocations(program_number, query_info,
//                                   filter_maskloc, lookup_segments);
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
// ```c
// /* Increment offset by 1 extra byte for the sentinel NULLB
//    between frames. */
// offset += length + 1;
// frame_offsets[context+1] = offset;
// ```
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2797
// ```c
// Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
// ...
// loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
// ...
// if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
//    Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
// ```
// prepare_lookup_query calls x_prepare_lookup_query with x_with_karlin = true: every context gets
// its Karlin block.
fn prepare_lookup_query(
    queries: &[Vec<QueryFrame>],
    ideal_params: KarlinParams,
    std_comp: &[f64; LOOKUP_ALPHABET_SIZE],
    word_length: usize,
    check_ideal: bool,
    matrix: ScoringMatrix,
    bounds: (i32, i32),
) -> PreparedLookupQuery {
    x_prepare_lookup_query(
        queries,
        ideal_params,
        std_comp,
        word_length,
        check_ideal,
        matrix,
        bounds,
        true,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2797
// ```c
// Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
// ...
// loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
// ...
// if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
//    Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
// ```
// NCBI computes the per-context Karlin blocks from the query composition, separately from the
// lookup table (the table is built in lookup_wrap.c from the query and its lookup segments only).
// This function runs the same preparation (prepare_lookup_query) and returns only the contexts.
/// The contexts `build_ncbi_lookup_for_profile` returns, for a caller that
/// drops the table (graduated 2026-10-22, was `LOSAT_X_LUTSPLIT`). The table
/// build only reads the contexts (their count and the end of the last one) and
/// hands them back unchanged, so they are those of `prepare_lookup_query`.
/// `None` when the build would print its diagnostics.
pub(crate) fn x_lookup_contexts_for_profile(
    queries: &[Vec<QueryFrame>],
    check_ideal: bool,
    matrix: ScoringMatrix,
    word_length: usize,
) -> Option<Vec<QueryContext>> {
    if diagnostics_enabled() {
        return None;
    }
    let bounds = if matrix == ScoringMatrix::Blosum62 {
        (-4, 11)
    } else {
        matrix_score_bounds(matrix)
    };
    Some(
        prepare_lookup_query(
            queries,
            ideal_karlin_params_for_matrix(matrix, bounds),
            &compute_std_aa_composition(),
            word_length,
            check_ideal,
            matrix,
            bounds,
        )
        .contexts,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803 (cited at
// `x_context_karlin`)
// NCBI computes the contexts' Karlin blocks one after the other; each reads only its own context's
// letters (`x_context_karlin`), so computing them on several threads gives each context the same
// block (as `x_prepare_blosum62_lookup_query_in_pool`). Scheduling only.
/// `x_lookup_contexts_for_profile` with the contexts' Karlin blocks computed on `pool`'s threads;
/// the same contexts.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) fn x_lookup_contexts_for_profile_in_pool(
    queries: &[Vec<QueryFrame>],
    check_ideal: bool,
    matrix: ScoringMatrix,
    word_length: usize,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Option<Vec<QueryContext>> {
    use rayon::prelude::*;
    if diagnostics_enabled() {
        return None;
    }
    let bounds = if matrix == ScoringMatrix::Blosum62 {
        (-4, 11)
    } else {
        matrix_score_bounds(matrix)
    };
    let ideal_params = ideal_karlin_params_for_matrix(matrix, bounds);
    let std_comp = compute_std_aa_composition();
    let mut contexts = x_prepare_lookup_query(
        queries,
        ideal_params,
        &std_comp,
        word_length,
        check_ideal,
        matrix,
        bounds,
        false,
    )
    .contexts;
    // The frames in the order of the contexts (`x_prepare_lookup_query`'s loops).
    let frames: Vec<&QueryFrame> = queries.iter().flatten().collect();
    let karlin: Vec<(bool, KarlinParams)> = pool.install(|| {
        frames
            .par_iter()
            .map(|frame| {
                x_context_karlin(frame, &std_comp, bounds, matrix, check_ideal, ideal_params)
            })
            .collect()
    });
    assert_eq!(karlin.len(), contexts.len());
    for (context, (is_valid, karlin_params)) in contexts.iter_mut().zip(karlin) {
        context.is_valid = is_valid;
        context.karlin_params = karlin_params;
    }
    Some(contexts)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
// ```c
// BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
//                       &lookup_wrap->lut);
// ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
// BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
//                           query, lookup_segments, 0);
// ```
// The table is built from the concatenated query, its lookup segments, the matrix and the
// threshold; it reads no Karlin parameter. This function skips the Karlin computation and builds
// the same table (the contexts it fills carry default Karlin parameters, which the table build
// does not read).
/// The table `build_ncbi_lookup_for_profile` returns, for a caller that drops
/// the contexts (graduated 2026-10-22, was `LOSAT_X_LUTSPLIT`). The
/// per-context Karlin parameters are the only part of the preparation the
/// table does not read; they are left out.
pub(crate) fn x_lookup_table_for_profile(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
) -> BlastAaLookupTable {
    x_lookup_table_for_profile_in(queries, threshold, matrix, word_length, None)
}

/// `x_lookup_table_for_profile`, the lean build's lists made on `pool`'s threads when it has two or
/// more (`x_lut_direct::x_build_lean_in`); the same table.
pub(crate) fn x_lookup_table_for_profile_in(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> BlastAaLookupTable {
    let bounds = if matrix == ScoringMatrix::Blosum62 {
        (-4, 11)
    } else {
        matrix_score_bounds(matrix)
    };
    let prepared = x_prepare_lookup_query(
        queries,
        ideal_karlin_params_for_matrix(matrix, bounds),
        &compute_std_aa_composition(),
        word_length,
        false,
        matrix,
        bounds,
        false,
    );
    // Dispatch point (LOSAT_X_TBNLUTLEAN, in the strict set): the same table without chains, every
    // temporary freed by the build (`x_lut_direct::x_build_lean`); in diagnostics mode, whose prints
    // read the chains, and with `=0`, the chain builder.
    let prepared = if x_tbn_lut_lean() {
        match x_lut_direct::x_build_lean_in(prepared, threshold, matrix, word_length, pool) {
            Ok(table) => {
                x_count_tbn_lut(0);
                return table;
            }
            Err(prepared) => prepared,
        }
    } else {
        prepared
    };
    x_count_tbn_lut(1);
    build_lookup_from_prepared(prepared, threshold, matrix, word_length).0
}

// No NCBI counterpart (memory only; the table is the one `build_lookup_from_prepared` builds):
// LOSAT_X_TBNLUTLEAN, in the strict set (on unless `LOSAT_X_TBNLUTLEAN=0`).
/// Whether TBLASTN builds its lookup tables with `x_lut_direct::x_build_lean`.
fn x_tbn_lut_lean() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNLUTLEAN"))
}

// No NCBI counterpart: LOSAT_X_STATS counters; they do not change any value NCBI computes.
/// TBLASTN tables built by `x_build_lean` (`[0]`) and by the chain builder (`[1]`).
static X_TBN_LUT_BUILDS: [std::sync::atomic::AtomicU64; 2] =
    [const { std::sync::atomic::AtomicU64::new(0) }; 2];

/// Counts one TBLASTN table built by `x_build_lean` (`which` 0) or the chain builder (1) for the
/// `[X_STATS] TBNLUTLEAN` line (`LOSAT_X_STATS` set).
fn x_count_tbn_lut(which: usize) {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if *ON.get_or_init(|| std::env::var_os("LOSAT_X_STATS").is_some()) {
        X_TBN_LUT_BUILDS[which].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }
}

/// The `[X_STATS] TBNLUTLEAN lean= chain=` line (printed by `main` at exit) with `LOSAT_X_STATS`,
/// when TBLASTN built any table.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn x_tbn_lut_print_summary() {
    let [lean, chain] = X_TBN_LUT_BUILDS
        .each_ref()
        .map(|c| c.load(std::sync::atomic::Ordering::Relaxed));
    if lean + chain > 0 {
        eprintln!("[X_STATS] TBNLUTLEAN lean={lean} chain={chain}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:104-108
// ```c
//    case eCompressedAaLookupTable:
//       BlastCompressedAaLookupTableNew(query, lookup_segments,
//                          (BlastCompressedAaLookupTable* *) &(lookup_wrap->lut),
//                          lookup_options, sbp);
// ```
/// The concatenated query and the lookup segments that `x_lookup_table_for_profile` builds its
/// table from (BLOSUM62, no Karlin blocks), for a table of another kind: TBLASTN's compressed
/// table (`tblastn/search_seed.rs`).
pub(crate) fn x_lookup_query_for_profile(
    queries: &[Vec<QueryFrame>],
    word_length: usize,
) -> (Vec<u8>, Vec<(i32, i32)>) {
    let bounds = (-4, 11);
    let prepared = x_prepare_lookup_query(
        queries,
        ideal_karlin_params_for_matrix(ScoringMatrix::Blosum62, bounds),
        &compute_std_aa_composition(),
        word_length,
        false,
        ScoringMatrix::Blosum62,
        bounds,
        false,
    );
    (prepared.concat_query, prepared.lookup_locations)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2803
// ```c
// Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
// ...
// loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
// if (loop_status) {
//     contexts[context].is_valid = FALSE;
//     ...
//     continue;
// }
// if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
//     Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
// ```
/// A context's validity and Karlin block, from its frame alone (the statements of
/// `x_prepare_lookup_query`'s loop).
fn x_context_karlin(
    frame: &QueryFrame,
    std_comp: &[f64; LOOKUP_ALPHABET_SIZE],
    bounds: (i32, i32),
    matrix: ScoringMatrix,
    check_ideal: bool,
    ideal_params: KarlinParams,
) -> (bool, KarlinParams) {
    let ctx_comp = compute_aa_composition(&frame.aa_seq, frame.aa_len);
    // NCBI c++/src/algo/blast/core/blast_stat.c:2778-2803:
    // score frequencies use this score block's matrix and lo/hi scores.
    let sfp =
        compute_score_freq_profile_for_matrix(&ctx_comp, std_comp, bounds.0, bounds.1, matrix);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2803
    // ```c
    // loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    // if (loop_status) {
    //     contexts[context].is_valid = FALSE;
    //     sbp->sfp[context] = Blast_ScoreFreqFree(sbp->sfp[context]);
    //     sbp->kbp_std[context] = Blast_KarlinBlkFree(sbp->kbp_std[context]);
    //     continue;
    // }
    // if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
    //     Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
    // ```
    match compute_karlin_params_ungapped(&sfp) {
        Ok(computed_params) => {
            let final_params = if check_ideal {
                apply_check_ideal(computed_params, ideal_params)
            } else {
                computed_params
            };
            (true, final_params)
        }
        Err(_) => {
            // NCBI leaves no valid per-context Karlin block behind on failure.
            // Rust stores a placeholder and requires later stages to honor
            // `is_valid`, matching `query_info->contexts[context].is_valid`.
            (false, KarlinParams::default())
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2797
// ```c
// Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
// ...
// loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
// ...
// if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
//    Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
// ```
// This is prepare_lookup_query with a flag: x_with_karlin = true runs the code of the reference
// path unchanged (the Karlin computation of `x_context_karlin`); false skips it (for
// `x_lookup_table_for_profile` and for `x_prepare_blosum62_lookup_query_in_pool`, which computes it
// on a pool).
fn x_prepare_lookup_query(
    queries: &[Vec<QueryFrame>],
    ideal_params: KarlinParams,
    std_comp: &[f64; LOOKUP_ALPHABET_SIZE],
    word_length: usize,
    check_ideal: bool,
    matrix: ScoringMatrix,
    bounds: (i32, i32),
    x_with_karlin: bool,
) -> PreparedLookupQuery {
    let mut prepared = PreparedLookupQuery {
        concat_query: Vec::new(),
        lookup_locations: Vec::new(),
        frame_bases: Vec::new(),
        contexts: Vec::new(),
        skipped_seg_mask: 0,
    };
    let mut query_offset = 0i32;

    for (q_idx, frames) in queries.iter().enumerate() {
        for (f_idx, frame) in frames.iter().enumerate() {
            prepared.frame_bases.push(query_offset);
            prepared.concat_query.extend_from_slice(&frame.aa_seq[1..]);

            let lookup_segments =
                build_lookup_segments_for_context(query_offset, frame.aa_len, &frame.seg_masks);
            let total_aa_positions = frame.aa_len.saturating_sub(word_length - 1);
            let lookup_positions = count_lookup_positions(&lookup_segments, word_length);
            prepared.skipped_seg_mask += total_aa_positions.saturating_sub(lookup_positions);
            prepared
                .lookup_locations
                .extend_from_slice(&lookup_segments);

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
            // ```c
            // BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
            //                       &lookup_wrap->lut);
            // ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
            // BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
            //                           query, lookup_segments, 0);
            // ```
            // The table build reads only the query, the lookup segments, the matrix and the threshold, so
            // the contexts can skip the Karlin parameters here.
            // A caller that only wants the table (`x_lookup_table_for_profile`).
            if !x_with_karlin {
                prepared.contexts.push(QueryContext {
                    q_idx: q_idx as u32,
                    f_idx: f_idx as u8,
                    frame: frame.frame,
                    aa_seq: frame.aa_seq.clone(),
                    aa_seq_nomask: frame.aa_seq_nomask.clone(),
                    aa_len: frame.aa_len,
                    orig_len: frame.orig_len,
                    frame_base: query_offset,
                    is_valid: true,
                    karlin_params: KarlinParams::default(),
                });
                query_offset += frame.aa_seq.len() as i32 - 1;
                continue;
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2803
            // ```c
            // loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
            // ```
            let (is_valid, final_params) =
                x_context_karlin(frame, std_comp, bounds, matrix, check_ideal, ideal_params);

            prepared.contexts.push(QueryContext {
                q_idx: q_idx as u32,
                f_idx: f_idx as u8,
                frame: frame.frame,
                aa_seq: frame.aa_seq.clone(),
                aa_seq_nomask: frame.aa_seq_nomask.clone(),
                aa_len: frame.aa_len,
                orig_len: frame.orig_len,
                frame_base: query_offset,
                is_valid,
                karlin_params: final_params,
            });

            query_offset += frame.aa_seq.len() as i32 - 1;
        }
    }

    prepared
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:621-636
// ```c
// if (lookup_segments) {
//     BLAST_ComplementMaskLocations(program_number, query_info,
//                                   filter_maskloc, lookup_segments);
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1349-1351
// ```c
// s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                 query, locations);
// ```
pub(crate) fn prepare_blosum62_lookup_query_for_word_size(
    queries: &[Vec<QueryFrame>],
    word_length: usize,
    check_ideal: bool,
) -> (Vec<u8>, Vec<(i32, i32)>, Vec<QueryContext>) {
    let ideal_params = blosum62_ideal_karlin_params();
    let std_comp = compute_std_aa_composition();
    let prepared = prepare_lookup_query(
        queries,
        ideal_params,
        &std_comp,
        word_length,
        check_ideal,
        ScoringMatrix::Blosum62,
        (-4, 11),
    );
    (
        prepared.concat_query,
        prepared.lookup_locations,
        prepared.contexts,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803
// ```c
//    for (context = query_info->first_context;
//         context <= query_info->last_context; ++context) {
//       ...
//       if ( !contexts[context].is_valid )
//           continue;
//       ...
//       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
// ```
// NCBI computes the contexts' Karlin blocks one after the other; each reads only its own context's
// letters (`x_context_karlin`), so computing them on several threads gives each context the same
// block. Scheduling only.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): `prepare_blosum62_lookup_query_for_word_size` with the contexts'
/// Karlin blocks (`x_context_karlin`) computed on `pool`'s threads; the same three values.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) fn x_prepare_blosum62_lookup_query_in_pool(
    queries: &[Vec<QueryFrame>],
    word_length: usize,
    check_ideal: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> (Vec<u8>, Vec<(i32, i32)>, Vec<QueryContext>) {
    use rayon::prelude::*;
    let lookup_query = XBlosum62LookupQuery::new(check_ideal);
    let (concat_query, lookup_locations, mut contexts) =
        lookup_query.without_karlin(queries, word_length);
    // The frames in the order of the contexts (`x_prepare_lookup_query`'s loops).
    let frames: Vec<&QueryFrame> = queries.iter().flatten().collect();
    let karlin: Vec<(bool, KarlinParams)> = pool.install(|| {
        frames
            .par_iter()
            .map(|frame| lookup_query.context_karlin(frame))
            .collect()
    });
    assert_eq!(karlin.len(), contexts.len());
    for (context, (is_valid, karlin_params)) in contexts.iter_mut().zip(karlin) {
        context.is_valid = is_valid;
        context.karlin_params = karlin_params;
    }
    (concat_query, lookup_locations, contexts)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803 (above)
// The values of `prepare_blosum62_lookup_query_for_word_size`, made in two steps: the lookup query,
// whose contexts carry no Karlin block yet, and each context's block (`x_context_karlin` with the
// arguments of `prepare_blosum62_lookup_query_for_word_size`), which reads only its own context's
// letters, so the caller may compute the blocks in any order and on any thread. Scheduling only.
/// EXPERIMENT (LOSAT_X_BPQBLOCK, LOSAT_X_BPLKPFLAT): the lookup query of
/// `prepare_blosum62_lookup_query_for_word_size` with the contexts' Karlin blocks left to the caller.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) struct XBlosum62LookupQuery {
    ideal_params: KarlinParams,
    std_comp: [f64; LOOKUP_ALPHABET_SIZE],
    check_ideal: bool,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl XBlosum62LookupQuery {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2833-2852 (the ideal block, as
    // `prepare_blosum62_lookup_query_for_word_size` reads it)
    // ```c
    //    Blast_ResFreqStdComp(sbp, stdrfp);
    //    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
    //    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
    // ```
    pub(crate) fn new(check_ideal: bool) -> Self {
        Self {
            ideal_params: blosum62_ideal_karlin_params(),
            std_comp: compute_std_aa_composition(),
            check_ideal,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:634-637 (the lookup segments,
    // `prepare_blosum62_lookup_query_for_word_size`)
    // ```c
    //     if (lookup_segments) {
    //         BLAST_ComplementMaskLocations(program_number, query_info,
    //                                       filter_maskloc, lookup_segments);
    //     }
    // ```
    /// The query, the lookup segments and the contexts of
    /// `prepare_blosum62_lookup_query_for_word_size`, each context valid with the default Karlin
    /// parameters until the caller sets them from `context_karlin` of its frame.
    pub(crate) fn without_karlin(
        &self,
        queries: &[Vec<QueryFrame>],
        word_length: usize,
    ) -> (Vec<u8>, Vec<(i32, i32)>, Vec<QueryContext>) {
        let prepared = x_prepare_lookup_query(
            queries,
            self.ideal_params,
            &self.std_comp,
            word_length,
            self.check_ideal,
            ScoringMatrix::Blosum62,
            (-4, 11),
            false,
        );
        (
            prepared.concat_query,
            prepared.lookup_locations,
            prepared.contexts,
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2803 (`x_context_karlin`)
    // ```c
    //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    // ```
    /// The validity and Karlin block of the context of `frame`.
    pub(crate) fn context_karlin(&self, frame: &QueryFrame) -> (bool, KarlinParams) {
        x_context_karlin(
            frame,
            &self.std_comp,
            (-4, 11),
            ScoringMatrix::Blosum62,
            self.check_ideal,
            self.ideal_params,
        )
    }
}

#[inline(always)]
fn lookup_chain_num_used(chain: &[i32]) -> usize {
    usize::try_from(chain[1]).expect("NCBI BLAST lookup chain hit count must fit in usize")
}

#[inline(always)]
fn lookup_chain_entries(chain: &[i32]) -> &[i32] {
    let num_used = lookup_chain_num_used(chain);
    &chain[2..2 + num_used]
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:33-77
// ```c
// index = ComputeTableIndex(wordsize, charsize, seq);
// if (backbone[index] == NULL) {
//     chain_size = 8;
//     hits_in_chain = 0;
//     chain = (Int4 *) malloc(chain_size * sizeof(Int4));
//     chain[0] = chain_size;
//     chain[1] = hits_in_chain;
//     backbone[index] = chain;
// }
// if ((hits_in_chain + 2) == chain_size) {
//     chain_size = chain_size * 2;
//     chain = (Int4 *) realloc(chain, chain_size * sizeof(Int4));
//     backbone[index] = chain;
//     chain[0] = chain_size;
// }
// chain[chain[1] + 2] = query_offset;
// chain[1]++;
// ```
#[inline]
fn blast_lookup_add_word_hit(
    backbone: &mut LookupBackboneChains,
    wordsize: usize,
    charsize: usize,
    seq: &[u8],
    query_offset: i32,
) {
    let mut index = 0usize;
    for &residue in seq.iter().take(wordsize) {
        index = (index << charsize) | residue as usize;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:48-72
    // ```c
    // if (backbone[index] == NULL) {
    //     chain_size = 8;
    //     hits_in_chain = 0;
    //     chain = (Int4 *) malloc(chain_size * sizeof(Int4));
    // ...
    // if ((hits_in_chain + 2) == chain_size) {
    //     chain_size = chain_size * 2;
    //     chain = (Int4 *) realloc(chain, chain_size * sizeof(Int4));
    // ...
    //     backbone[index] = chain;
    //     chain[0] = chain_size;
    // ```
    // The chain is created and grown through x_zeroed8 and x_resize, which use the arena. The
    // statements that follow are the C ones.
    let LookupBackboneChains { cells, arena } = backbone;
    let chain = cells[index].get_or_insert_with(|| {
        let mut chain = LookupChain::x_zeroed8(arena);
        chain[0] = 8;
        chain[1] = 0;
        chain
    });

    let hits_in_chain = lookup_chain_num_used(chain);
    let chain_size = usize::try_from(chain[0]).expect("NCBI BLAST lookup chain size must fit");
    if (hits_in_chain + 2) == chain_size {
        let new_size = chain_size * 2;
        chain.x_resize(new_size, arena);
        chain[0] = i32::try_from(new_size).expect("NCBI BLAST lookup chain size must fit Int4");
    }

    chain[hits_in_chain + 2] = query_offset;
    chain[1] += 1;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:79-132
// ```c
// for (loc = locations; loc; loc = loc->next) {
//     Int4 from = loc->ssr->left;
//     Int4 to = loc->ssr->right;
//     if (word_length > to - from + 1)
//         continue;
//     seq = query->sequence + from;
//     word_target = seq + lut_word_length;
//     for (offset = from; offset <= to; offset++, seq++) {
//         if (seq >= word_target) {
//             BlastLookupAddWordHit(backbone, lut_word_length, charsize,
//                                   seq - lut_word_length, offset - lut_word_length);
//         }
//         if (*seq & invalid_mask)
//             word_target = seq + lut_word_length + 1;
//     }
//     if (seq >= word_target) {
//         BlastLookupAddWordHit(backbone, lut_word_length, charsize,
//                               seq - lut_word_length, offset - lut_word_length);
//     }
// }
// ```
fn blast_lookup_index_query_exact_matches(
    backbone: &mut LookupBackboneChains,
    word_length: i32,
    charsize: i32,
    lut_word_length: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    skipped_invalid_residue: &mut usize,
) -> usize {
    let invalid_mask: u8 = 0xffu8 << (charsize as u32);
    let mut total_exact_positions = 0usize;

    for &(from, to) in locations {
        if word_length > to - from + 1 {
            continue;
        }

        let mut offset = from;
        let mut seq = usize::try_from(from).expect("NCBI BLAST query offset must fit usize");
        let mut word_target = seq + usize::try_from(lut_word_length).expect("word length fits");

        while offset <= to {
            if seq >= word_target {
                let word_start = seq - usize::try_from(lut_word_length).expect("word length fits");
                blast_lookup_add_word_hit(
                    backbone,
                    usize::try_from(lut_word_length).expect("word length fits"),
                    usize::try_from(charsize).expect("charsize fits"),
                    &query[word_start..],
                    offset - lut_word_length,
                );
                total_exact_positions += 1;
            }

            if (query[seq] & invalid_mask) != 0 {
                word_target = seq + usize::try_from(lut_word_length).expect("word length fits") + 1;
                *skipped_invalid_residue += 1;
            }
            offset += 1;
            seq += 1;
        }

        if seq >= word_target {
            let word_start = seq - usize::try_from(lut_word_length).expect("word length fits");
            blast_lookup_add_word_hit(
                backbone,
                usize::try_from(lut_word_length).expect("word length fits"),
                usize::try_from(charsize).expect("charsize fits"),
                &query[word_start..],
                offset - lut_word_length,
            );
            total_exact_positions += 1;
        }
    }

    total_exact_positions
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:546-606
// ```c
// score -= info->row_max[query_word[current_pos]];
// row = info->matrix[query_word[current_pos]];
// if (current_pos == info->wordsize - 1) {
//     for (i = 0; i < alphabet_size; i++) {
//         if (score + row[i] >= threshold) {
//             subject_word[current_pos] = i;
//             for (j = 0; j < offset_list[1]; j++) {
//                 BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
//                                       charsize, subject_word,
//                                       query_bias + offset_list[j + 2]);
//             }
//         }
//     }
//     return;
// }
// for (i = 0; i < alphabet_size; i++) {
//     if (score + row[i] >= threshold) {
//         subject_word[current_pos] = i;
//         s_AddWordHitsCore(info, score + row[i], current_pos + 1);
//     }
// }
// ```
fn blast_aa_add_word_hits_core(info: &mut NeighborInfo<'_>, score: i32, current_pos: usize) {
    let query_residue = info.query_word[current_pos] as usize;
    let score = score - info.row_max[query_residue];
    // NCBI c++/src/algo/blast/core/blast_aalookup.c:562-563:
    // row = info->matrix[query_word[current_pos]];
    let matrix = info.matrix;

    if current_pos == info.wordsize - 1 {
        for residue in 0..info.alphabet_size {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:580-582
            // ```c
            // for (i = 0; i < alphabet_size; i++) {
            //     if (score + row[i] >= threshold) {
            // ```
            let residue_score = lookup_matrix_score(matrix, query_residue, residue);
            if score + residue_score >= info.threshold {
                info.subject_word[current_pos] = residue as u8;
                for &query_offset in lookup_chain_entries(info.offset_list) {
                    blast_lookup_add_word_hit(
                        info.thin_backbone,
                        info.wordsize,
                        info.charsize,
                        &info.subject_word,
                        info.query_bias + query_offset,
                    );
                    info.neighbor_added_count += 1;
                }
                info.neighbor_words_generated += 1;
            }
        }
        return;
    }

    for residue in 0..info.alphabet_size {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:600-603
        // ```c
        // for (i = 0; i < alphabet_size; i++) {
        //     if (score + row[i] >= threshold) {
        //         subject_word[current_pos] = i;
        //         s_AddWordHitsCore(info, score + row[i], current_pos + 1);
        // ```
        let residue_score = lookup_matrix_score(matrix, query_residue, residue);
        if score + residue_score >= info.threshold {
            info.subject_word[current_pos] = residue as u8;
            blast_aa_add_word_hits_core(info, score + residue_score, current_pos + 1);
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:472-544
// ```c
// w = query + offset_list[2];
// score = matrix[w[0]][w[0]];
// for (i = 1; i < lookup->word_length; i++)
//     score += matrix[w[i]][w[i]];
// if (lookup->threshold == 0 || score < lookup->threshold) {
//     for (i = 0; i < offset_list[1]; i++) {
//         BlastLookupAddWordHit(lookup->thin_backbone, lookup->word_length,
//                               lookup->charsize, w,
//                               query_bias + offset_list[i + 2]);
//     }
// }
// if (lookup->threshold == 0)
//     return;
// score = row_max[w[0]];
// for (i = 1; i < lookup->word_length; i++)
//     score += row_max[w[i]];
// s_AddWordHitsCore(&info, score, 0);
// ```
fn blast_aa_add_word_hits(
    thin_backbone: &mut LookupBackboneChains,
    alphabet_size: usize,
    word_length: usize,
    charsize: usize,
    threshold: i32,
    query: &[u8],
    offset_list: &[i32],
    query_bias: i32,
    row_max: &[i32],
    matrix: ScoringMatrix,
) -> (usize, usize, usize) {
    let query_offset = usize::try_from(offset_list[2]).expect("NCBI BLAST query offset must fit");
    let query_word = &query[query_offset..];
    let query_residue = query_word[0] as usize;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:492-496
    // ```c
    // score = matrix[w[0]][w[0]];
    // for (i = 1; i < lookup->word_length; i++)
    //     score += matrix[w[i]][w[i]];
    // ```
    // NCBI c++/src/algo/blast/core/blast_aalookup.c:492-496:
    // score = matrix[w[0]][w[0]]; then add each matrix[w[i]][w[i]].
    let mut self_score = lookup_matrix_score(matrix, query_residue, query_residue);
    for residue in query_word.iter().take(word_length).skip(1) {
        let residue = *residue as usize;
        self_score += lookup_matrix_score(matrix, residue, residue);
    }

    let mut exact_added_count = 0usize;
    if threshold == 0 || self_score < threshold {
        for &query_offset in lookup_chain_entries(offset_list) {
            blast_lookup_add_word_hit(
                thin_backbone,
                word_length,
                charsize,
                query_word,
                query_bias + query_offset,
            );
            exact_added_count += 1;
        }
    }
    if threshold == 0 {
        return (exact_added_count, 0, 0);
    }

    let mut info = NeighborInfo {
        thin_backbone,
        query_word,
        subject_word: [0; 32],
        alphabet_size,
        wordsize: word_length,
        charsize,
        row_max,
        offset_list,
        threshold,
        query_bias,
        matrix,
        neighbor_added_count: 0,
        neighbor_words_generated: 0,
    };
    let mut max_score = row_max[query_word[0] as usize];
    for residue in query_word.iter().take(word_length).skip(1) {
        max_score += row_max[*residue as usize];
    }
    blast_aa_add_word_hits_core(&mut info, max_score, 0);

    (
        exact_added_count,
        info.neighbor_added_count,
        info.neighbor_words_generated,
    )
}

/// Build NCBI-style lookup table with presence vector and overflow.
///
/// Reference: blast_aalookup.c BlastAaLookupIndexQuery + BlastAaLookupFinalize
///
/// Neighbor words are always precomputed in the lookup table (no lazy mode).
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-543
///
/// NCBI reference (fixed alphabet and neighbor precomputation):
/// ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:436-456
/// ```c
/// ASSERT(lookup->alphabet_size <= BLASTAA_SIZE);
/// ...
/// BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
///                                   lookup->charsize, lookup->word_length,
///                                   query, location);
/// ```
pub fn build_ncbi_lookup(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    karlin_params: &crate::stats::KarlinParams,
    check_ideal: bool,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    build_ncbi_lookup_for_profile(
        queries,
        threshold,
        karlin_params,
        check_ideal,
        ScoringMatrix::Blosum62,
        LOOKUP_WORD_LENGTH,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
// ```c
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (exact_backbone[i] != NULL) {
//             s_AddWordHits(lookup, matrix, query->sequence,
//                           exact_backbone[i], query_bias, row_max);
//             sfree(exact_backbone[i]);
//         }
//     }
// ```
/// [`build_ncbi_lookup`] for BLASTP: the table is filled in place from the per-cell counts, in
/// the same append order as the loop above and BlastAaLookupFinalize (x_lut_direct.rs), so the
/// chains of the reference builder are never held.
pub(crate) fn x_bp_build_ncbi_lookup(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    check_ideal: bool,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    x_lut_direct::x_bp_build_lookup(queries, threshold, check_ideal)
}

/// EXPERIMENT (LOSAT_X_BPLK16): [`x_bp_build_ncbi_lookup`] with the width of the table's hits
/// chosen by the build (`narrow`: a table of `u16` hits where the lookup query is at most 65,535
/// letters and the in-place builder applies; x_lut_direct.rs `x_bp_build_lookup_as`).
pub(crate) fn x_bp_build_ncbi_lookup_as(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    check_ideal: bool,
    narrow: bool,
) -> (x_lut_direct::XBpLookup, Vec<QueryContext>) {
    x_lut_direct::x_bp_build_lookup_as(queries, threshold, check_ideal, narrow)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:356-360
// ```c
// for (j=0; j <lookup->thin_backbone[i][1]; j++)
//     dest[j] = lookup->thin_backbone[i][j + 2];
// /* done with this chain- free it */
// sfree(lookup->thin_backbone[i]);
// lookup->thin_backbone[i] = NULL;
// ```
/// EXPERIMENT (LOSAT_X_TBXLUT): [`build_ncbi_lookup`] for the full table
/// of a TBLASTX query batch. With the switch the table is filled in place from the per-cell counts,
/// in the same append order as the reference builder and BlastAaLookupFinalize (x_lut_direct.rs),
/// so no chain is held: NCBI frees each chain as it copies it (above), the reference builder here
/// holds all of them, and the arena's blocks, until the table is complete.
pub(crate) fn x_tbx_build_ncbi_lookup(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    karlin_params: &crate::stats::KarlinParams,
    check_ideal: bool,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    x_lut_direct::x_tbx_build_lookup(queries, threshold, karlin_params, check_ideal)
}

// NCBI c++/src/algo/blast/core/blast_aalookup.c:227-245,438-465:
// lookup->word_length = opt->word_size; lookup->threshold = opt->threshold;
// row_max[i] = MAX(row_max[i], matrix[i][j]);
// s_AddWordHits(lookup, matrix, query->sequence, ...);
#[allow(dead_code)]
pub(crate) fn build_ncbi_lookup_for_profile(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    _karlin_params: &crate::stats::KarlinParams,
    check_ideal: bool,
    matrix: ScoringMatrix,
    word_length: usize,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    let bounds = if matrix == ScoringMatrix::Blosum62 {
        (-4, 11)
    } else {
        matrix_score_bounds(matrix)
    };
    let prepared = prepare_lookup_query(
        queries,
        ideal_karlin_params_for_matrix(matrix, bounds),
        &compute_std_aa_composition(),
        word_length,
        check_ideal,
        matrix,
        bounds,
    );
    build_lookup_from_prepared(prepared, threshold, matrix, word_length)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469
// ```c++
//     /* create an empty backbone */
//
//     exact_backbone = (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
//
//     /* find all the exact matches, grouping together all offsets of identical
//        query words. The query bias is not used here, since the next stage
//        will need real offsets into the query sequence */
//
//     BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
//                                       lookup->charsize, lookup->word_length,
//                                       query, location);
//
//     /* walk though the list of exact matches previously computed. Find
//        neighboring words for entire lists at a time */
//
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (exact_backbone[i] != NULL) {
//             s_AddWordHits(lookup, matrix, query->sequence,
//                           exact_backbone[i], query_bias, row_max);
//             sfree(exact_backbone[i]);
//         }
//     }
//
//     sfree(exact_backbone);
// ```
/// BLASTX's lookup table of a prepared query; `parts` above 1: a pooled build's passes on that many tasks of the pool it runs on
/// (`LOSAT_X_BXBUILDPAR`, `x_lut_direct.rs` `x_partitioned_passes`); the table is the same.
pub(crate) fn build_ncbi_lookup_from_prepared(
    concat_query: Vec<u8>,
    lookup_locations: Vec<(i32, i32)>,
    contexts: Vec<QueryContext>,
    threshold: i32,
    pooled: bool,
    parts: usize,
) -> BlastAaLookupTable {
    let frame_bases = contexts.iter().map(|c| c.frame_base).collect();
    let prepared = PreparedLookupQuery {
        concat_query,
        lookup_locations,
        contexts,
        frame_bases,
        skipped_seg_mask: 0,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
    // ```c
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    //             sfree(exact_backbone[i]);
    //         }
    //     }
    // ```
    // BLASTX: the table is filled in place from the per-cell counts, in the same append order as
    // this loop and BlastAaLookupFinalize (x_lut_direct.rs); where that declines (diagnostics
    // mode, or a backbone index wider than 16 bits) the reference builder below runs. With `pooled`
    // (the pooled build, chosen by the caller for its search, LOSAT_X_BXPOOLN1) the build takes
    // kept vectors; the caller gives them back (`x_give_back_pooled`) when it drops the table, or
    // the pool goes without them (a missed reuse only).
    let built = if pooled {
        x_lut_direct::x_build_direct_pooled(prepared, threshold, ScoringMatrix::Blosum62, 3, parts)
    } else {
        x_lut_direct::x_build_direct(prepared, threshold, ScoringMatrix::Blosum62, 3)
    };
    let prepared = match built {
        Ok(table) => return table,
        Err(prepared) => prepared,
    };
    build_lookup_from_prepared(prepared, threshold, ScoringMatrix::Blosum62, 3).0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:227-245
// ```c++
// Int4 BlastAaLookupTableNew(const LookupTableOptions * opt,
//                            BlastAaLookupTable * *lut)
// {
//     Int4 i;
//     BlastAaLookupTable *lookup = *lut =
//         (BlastAaLookupTable *) calloc(1, sizeof(BlastAaLookupTable));
//
//     ASSERT(lookup != NULL);
//
//     lookup->charsize = ilog2(BLASTAA_SIZE) + 1;
//     lookup->word_length = opt->word_size;
//
//     for (i = 0; i < lookup->word_length; i++)
//         lookup->backbone_size |= (BLASTAA_SIZE - 1) << (i * lookup->charsize);
//     lookup->backbone_size++;
//
//     lookup->mask = (1 << (opt->word_size * lookup->charsize)) - 1;
//     lookup->alphabet_size = BLASTAA_SIZE;
//     lookup->threshold = (Int4)opt->threshold;
// ```
fn build_lookup_from_prepared(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    let diag_enabled = diagnostics_enabled();
    let alphabet_size = LOOKUP_ALPHABET_SIZE;
    let charsize = ilog2(alphabet_size) + 1;
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);
    let build_stats_initial = prepared.skipped_seg_mask;
    let mut build_stats = LookupBuildStats {
        skipped_seg_mask: build_stats_initial,
        ..LookupBuildStats::default()
    };
    let PreparedLookupQuery {
        concat_query,
        lookup_locations,
        frame_bases,
        contexts,
        skipped_seg_mask: _,
    } = prepared;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:246-250
    // ```c
    // Uint4 QueryInfo_GetSeqBufLen(const BlastQueryInfo* qinfo)
    // {
    //     BlastContextInfo * cinfo = & qinfo->contexts[qinfo->last_context];
    //     return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:500-501,654
    // ```c
    // int buflen = QueryInfo_GetSeqBufLen(qinfo);
    // TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
    // ...
    // BlastSeqBlkSetSequence(*seqblk, buf.release(), buflen - 2);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:147-160
    // ```c
    // seq_blk->sequence = (Uint1*) sequence + 1;
    // seq_blk->length = seqlen;
    // ```
    //
    // `concat_query` keeps the trailing NULLB so the scan/extension code can
    // read sentinel bytes just like NCBI's allocated buffer. However
    // `BLAST_SequenceBlk->length` excludes the final trailing NULLB, so
    // `lookup.query_length` must stop at the last real residue of the last
    // context, not at `concat_query.len()`.
    let query_length = contexts
        .last()
        .map(|ctx| ctx.frame_base + ctx.aa_len as i32)
        .unwrap_or(0);

    // NCBI c++/src/algo/blast/core/blast_aalookup.c:438-444:
    // row_max[i] = matrix[i][0];
    // for (j = 1; j < alphabet_size; j++) row_max[i] = MAX(row_max[i], matrix[i][j]);
    let row_max: Vec<i32> = (0..alphabet_size)
        .map(|i| {
            (0..alphabet_size)
                .map(|j| lookup_matrix_score(matrix, i, j))
                .max()
                .unwrap_or(i32::MIN)
        })
        .collect();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-469
    // ```c
    // exact_backbone = (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
    // BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
    //                                   lookup->charsize, lookup->word_length,
    //                                   query, location);
    // for (i = 0; i < lookup->backbone_size; i++) {
    //     if (exact_backbone[i] != NULL) {
    //         s_AddWordHits(lookup, matrix, query->sequence,
    //                       exact_backbone[i], query_bias, row_max);
    //     }
    // }
    // ```
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:448
    // ```c
    // exact_backbone = (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
    // ```
    // Same empty backbone as the C code; its chains come from an arena.
    let mut exact_backbone = LookupBackboneChains::new(backbone_size);
    build_stats.total_exact_positions = blast_lookup_index_query_exact_matches(
        &mut exact_backbone,
        word_length as i32,
        charsize as i32,
        word_length as i32,
        &concat_query,
        &lookup_locations,
        &mut build_stats.skipped_invalid_residue,
    );

    // Phase 1a diagnostics (guarded)
    if diag_enabled {
        let unique_exact_words = exact_backbone
            .iter()
            .filter(|chain| chain.is_some())
            .count();
        let max_exact_per_word = exact_backbone
            .iter()
            .filter_map(|chain| chain.as_deref())
            .map(lookup_chain_num_used)
            .max()
            .unwrap_or(0);
        eprintln!("\n=== Phase 1a: Exact Indexing Diagnostics ===");
        eprintln!(
            "Total exact positions indexed: {}",
            build_stats.total_exact_positions
        );
        eprintln!("Unique exact words: {}", unique_exact_words);
        eprintln!("Max offsets per exact word: {}", max_exact_per_word);
        eprintln!(
            "Skipped (invalid residue): {}",
            build_stats.skipped_invalid_residue
        );
        eprintln!("Skipped (SEG mask): {}", build_stats.skipped_seg_mask);
    }

    // Phase 1b: Add neighboring words (NCBI precomputed neighbors).
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-543
    let residue_mask: usize = (1usize << charsize) - 1;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:246-247
    // ```c
    // lookup->thin_backbone =
    //     (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
    // ```
    // Same empty thin backbone as the C code; its chains come from an arena.
    let mut thin_backbone = LookupBackboneChains::new(backbone_size);
    let mut counts: Vec<u32> = vec![0; backbone_size];

    for idx in 0..backbone_size {
        let Some(offset_list) = exact_backbone[idx].as_deref() else {
            continue;
        };
        let (exact_added_count, neighbor_added_count, neighbor_words_generated) =
            blast_aa_add_word_hits(
                &mut thin_backbone,
                alphabet_size,
                word_length,
                charsize,
                threshold,
                &concat_query,
                offset_list,
                0,
                &row_max,
                matrix,
            );
        build_stats.exact_added_count += exact_added_count;
        build_stats.neighbor_added_count += neighbor_added_count;
        build_stats.neighbor_words_generated += neighbor_words_generated;
        if threshold == 0 {
            build_stats.words_with_exact_only += 1;
        }
        if neighbor_words_generated > 0 {
            build_stats.words_with_neighbors += 1;
            if neighbor_words_generated > build_stats.max_neighbors_for_single_word {
                build_stats.max_neighbors_for_single_word = neighbor_words_generated;
                build_stats.max_neighbor_word_idx = idx;
            }
        }
    }
    for idx in 0..backbone_size {
        counts[idx] = thin_backbone[idx]
            .as_deref()
            .map(lookup_chain_num_used)
            .unwrap_or(0) as u32;
    }

    // Phase 1b diagnostics (guarded)
    if diag_enabled {
        eprintln!("\n=== Phase 1b: Neighbor Generation Diagnostics ===");
        eprintln!("Threshold: {}", threshold);
        eprintln!(
            "Exact entries added (self_score < threshold): {}",
            build_stats.exact_added_count
        );
        eprintln!(
            "Neighbor entries added: {}",
            build_stats.neighbor_added_count
        );
        eprintln!(
            "Neighbor words generated (unique): {}",
            build_stats.neighbor_words_generated
        );
        eprintln!(
            "Words with exact-only additions: {}",
            build_stats.words_with_exact_only
        );
        eprintln!(
            "Words with neighbor generation: {}",
            build_stats.words_with_neighbors
        );
        eprintln!(
            "Max neighbors generated for single word: {}",
            build_stats.max_neighbors_for_single_word
        );
        if build_stats.max_neighbors_for_single_word > 0 {
            let w0 = (build_stats.max_neighbor_word_idx >> (2 * charsize)) & residue_mask;
            let w1 = (build_stats.max_neighbor_word_idx >> charsize) & residue_mask;
            let w2 = build_stats.max_neighbor_word_idx & residue_mask;
            eprintln!(
                "  Word with max neighbors: idx={} (residues {},{},{})",
                build_stats.max_neighbor_word_idx, w0, w1, w2
            );
        }
        let expansion_factor = if build_stats.total_exact_positions > 0 {
            (build_stats.exact_added_count + build_stats.neighbor_added_count) as f64
                / build_stats.total_exact_positions as f64
        } else {
            0.0
        };
        eprintln!(
            "Expansion factor (total_entries / exact_positions): {:.2}x",
            expansion_factor
        );
    }

    // Phase 2: Finalize
    let mut backbone: Vec<BackboneCell> = vec![BackboneCell::default(); backbone_size];
    // NCBI allocates pv as: (backbone_size >> PV_ARRAY_BTS) + 1.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:327-329
    let pv_size = (backbone_size >> PV_ARRAY_BTS) + 1;
    let mut pv: Vec<u32> = vec![0u32; pv_size];

    // Pre-finalize diagnostics: analyze high-frequency k-mers (guarded)
    if diag_enabled {
        eprintln!("\n=== Phase 2: High-Frequency K-mer Analysis ===");

        // Collect (count, idx) pairs and sort by count descending
        let mut count_idx_pairs: Vec<(usize, usize)> = counts
            .iter()
            .enumerate()
            .filter(|(_, &c)| c > 0)
            .map(|(i, &c)| (c as usize, i))
            .collect();
        count_idx_pairs.sort_by(|a, b| b.0.cmp(&a.0));

        // Show top 10 high-frequency k-mers
        eprintln!("Top 10 high-frequency k-mers (before suppression):");
        for (rank, &(count, idx)) in count_idx_pairs.iter().take(10).enumerate() {
            let w0 = (idx >> (2 * charsize)) & residue_mask;
            let w1 = (idx >> charsize) & residue_mask;
            let w2 = idx & residue_mask;
            eprintln!(
                "  #{}: count={}, idx={} (residues {},{},{})",
                rank + 1,
                count,
                idx,
                w0,
                w1,
                w2
            );
        }

        // Distribution histogram
        let hist_buckets = [
            1,
            10,
            100,
            500,
            1000,
            5000,
            10000,
            50000,
            100000,
            usize::MAX,
        ];
        let mut hist: Vec<usize> = vec![0; hist_buckets.len()];
        let mut hist_entries: Vec<usize> = vec![0; hist_buckets.len()];
        for &(count, _) in &count_idx_pairs {
            for (i, &bucket) in hist_buckets.iter().enumerate() {
                if count <= bucket {
                    hist[i] += 1;
                    hist_entries[i] += count;
                    break;
                }
            }
        }
        eprintln!("Hit count distribution:");
        let mut prev = 0;
        for (i, &bucket) in hist_buckets.iter().enumerate() {
            if hist[i] > 0 {
                let label = if bucket == usize::MAX {
                    format!(">{}", prev)
                } else {
                    format!("{}-{}", prev + 1, bucket)
                };
                eprintln!(
                    "  {}: {} cells, {} entries",
                    label, hist[i], hist_entries[i]
                );
            }
            prev = bucket;
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:282-299
    // ```c
    // for (i = 0; i < lookup->backbone_size; i++) {
    //     if (lookup->thin_backbone[i]) {
    //         if (lookup->thin_backbone[i][1] > AA_HITS_PER_CELL)
    //             overflow_cells_needed += lookup->thin_backbone[i][1];
    //         if (lookup->thin_backbone[i][1] > longest_chain)
    //             longest_chain = lookup->thin_backbone[i][1];
    //     }
    // }
    // ```
    let mut overflow_size = 0usize;
    for &count in &counts {
        let count = count as usize;
        if count > AA_HITS_PER_CELL {
            overflow_size += count;
        }
    }

    let mut overflow: Vec<i32> = vec![0; overflow_size];
    let mut overflow_cursor = 0usize;

    let mut longest_chain: i32 = 0;
    for &count in &counts {
        if (count as i32) > longest_chain {
            longest_chain = count as i32;
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:337-360
    // ```c
    // for (i = 0; i < lookup->backbone_size; i++) {
    //     if (lookup->thin_backbone[i]) {
    //         PV_SET(pv, i, PV_ARRAY_BTS);
    //         bbc[i].num_used = lookup->thin_backbone[i][1];
    //         if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
    //             dest = bbc[i].payload.entries;
    //         else {
    //             bbc[i].payload.overflow_cursor = overflow_cursor;
    //             dest = (Int4 *) lookup->overflow + overflow_cursor;
    //             overflow_cursor += lookup->thin_backbone[i][1];
    //         }
    //         for (j = 0; j < lookup->thin_backbone[i][1]; j++)
    //             dest[j] = lookup->thin_backbone[i][j + 2];
    //     }
    // }
    // ```
    for idx in 0..backbone_size {
        let Some(chain) = thin_backbone[idx].as_deref() else {
            continue;
        };
        let entries = lookup_chain_entries(chain);
        let count = entries.len();
        pv_set(&mut pv, idx);
        backbone[idx].num_used = count as i32;

        if count <= AA_HITS_PER_CELL {
            for (i, &off) in entries.iter().enumerate() {
                backbone[idx].entries[i] = off;
            }
        } else {
            backbone[idx].entries[0] = overflow_cursor as i32;
            for &off in entries {
                overflow[overflow_cursor] = off;
                overflow_cursor += 1;
            }
        }
    }

    let nonempty = backbone.iter().filter(|c| c.num_used > 0).count();
    let total: usize = backbone.iter().map(|c| c.num_used.max(0) as usize).sum();
    if diag_enabled {
        eprintln!(
            "\n=== NCBI Backbone Lookup Table (BLASTAA_SIZE={}) ===",
            LOOKUP_ALPHABET_SIZE
        );
        eprintln!("Contexts: {}", contexts.len());
        eprintln!(
            "Backbone size: {} cells ({:.1} MB)",
            backbone_size,
            (backbone_size * std::mem::size_of::<BackboneCell>()) as f64 / 1e6
        );
        eprintln!("Non-empty cells: {}", nonempty);
        eprintln!("Total entries: {}", total);
        eprintln!(
            "Overflow size: {} ({:.1} MB)",
            overflow_cursor,
            (overflow_cursor * 4) as f64 / 1e6
        );
        eprintln!("Longest chain: {}", longest_chain);
        eprintln!("==========================================\n");
    }

    (
        BlastAaLookupTable {
            backbone,
            overflow,
            pv,
            frame_bases,
            num_contexts: contexts.len(),
            query_length,
            word_length: word_length as i32,
            alphabet_size: alphabet_size as i32,
            charsize: charsize as i32,
            mask: mask as i32,
            longest_chain,
            threshold,
            row_max,
        },
        contexts,
    )
}

/// Build a simple direct (exact) 3-mer lookup table for tests and diagnostics.
pub fn build_direct_lookup(queries: &[Vec<QueryFrame>]) -> Vec<Vec<(u32, u8, u32)>> {
    let word_length = LOOKUP_WORD_LENGTH;
    let alphabet_size = LOOKUP_ALPHABET_SIZE;
    let charsize = ilog2(alphabet_size) + 1;
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);

    let mut table: Vec<Vec<(u32, u8, u32)>> = vec![Vec::new(); backbone_size];

    for (q_idx, frames) in queries.iter().enumerate() {
        for (f_idx, frame) in frames.iter().enumerate() {
            if frame.aa_len < 3 || frame.aa_seq.len() < 5 {
                continue;
            }

            let seq_ptr = &frame.aa_seq[1..frame.aa_seq.len() - 1];
            // NCBI invalid_mask = 0xff << charsize (Uint1).
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:89-120
            let invalid_mask: u8 = 0xffu8 << (charsize as u32);
            let unmasked = compute_unmasked_intervals(&frame.seg_masks, frame.aa_len);
            for (start, end) in unmasked {
                let from = start as i32;
                let to = end.saturating_sub(1) as i32;
                // NCBI: if (word_length > to - from + 1) continue;
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:96-100
                if word_length as i32 > to - from + 1 {
                    continue;
                }
                let mut offset = from;
                let mut word_target = from + word_length as i32;

                while offset <= to {
                    if offset >= word_target {
                        let word_start = (offset - word_length as i32) as usize;
                        let idx = encode_kmer_3(
                            seq_ptr[word_start] as usize,
                            seq_ptr[word_start + 1] as usize,
                            seq_ptr[word_start + 2] as usize,
                            charsize,
                        ) & mask;
                        table[idx].push((q_idx as u32, f_idx as u8, word_start as u32));
                    }
                    if (seq_ptr[offset as usize] & invalid_mask) != 0 {
                        // NCBI: word_target = seq + lut_word_length + 1
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:117-120
                        word_target = offset + word_length as i32 + 1;
                    }
                    offset += 1;
                }

                if offset >= word_target {
                    let word_start = (offset - word_length as i32) as usize;
                    let idx = encode_kmer_3(
                        seq_ptr[word_start] as usize,
                        seq_ptr[word_start + 1] as usize,
                        seq_ptr[word_start + 2] as usize,
                        charsize,
                    ) & mask;
                    table[idx].push((q_idx as u32, f_idx as u8, word_start as u32));
                }
            }
        }
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ScoringMatrix;
    use crate::stats::lookup_protein_params_ungapped;
    use crate::utils::matrix::{blosum62_score, ncbistdaa, BLASTAA_SIZE};

    // No NCBI counterpart (LOSAT_X_BPDBKEEP): the unthresholded arena cuts its slabs and long chains
    // from one region while it has room (one region per 32 MiB on 64-bit glibc, else one per
    // allocation): the blocks it gives do not overlap and keep what was written to them.
    #[test]
    fn unthresholded_slabs_are_cut_from_one_region() {
        let mut arena = XChainArena::new();
        arena.x_unthresholded = true;
        let sizes = [
            XChainArena::SLAB_INTS,
            1 << 15,
            XChainArena::SLAB_INTS,
            3,
            1 << 17,
        ];
        let mut blocks = Vec::new();
        for (i, &ints) in sizes.iter().enumerate() {
            let start = arena.allocate(ints);
            // SAFETY: the arena gave `ints` ints at `start`.
            unsafe { std::slice::from_raw_parts_mut(start, ints) }.fill(i as i32 + 1);
            blocks.push((start as usize, ints));
        }
        for (i, &(start, ints)) in blocks.iter().enumerate() {
            // SAFETY: as above; the arena is alive.
            let ints_read = unsafe { std::slice::from_raw_parts(start as *const i32, ints) };
            assert!(ints_read.iter().all(|&v| v == i as i32 + 1), "block {i}");
            for &(other, other_ints) in &blocks[i + 1..] {
                let end = start + ints * 4;
                assert!(
                    end <= other || other + other_ints * 4 <= start,
                    "blocks {i} overlap"
                );
            }
        }
        let regions = if crate::utils::x_mapped::X_MAPPED_ALLOCATOR {
            1
        } else {
            sizes.len()
        };
        assert_eq!(arena.x_unthresholded_slabs.len(), regions);
        arena.x_reset();
        assert!(arena.x_unthresholded_slabs.is_empty());
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:87-129
    // ```c
    // if (seq >= word_target) {
    //     BlastLookupAddWordHit(backbone, lut_word_length, charsize,
    //                           seq - lut_word_length, offset - lut_word_length);
    // }
    // ```
    fn make_poly_a_query_frame() -> QueryFrame {
        QueryFrame {
            frame: 1,
            aa_seq: vec![
                0,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                0,
            ],
            aa_seq_nomask: None,
            aa_len: 6,
            orig_len: 6,
            seg_masks: Vec::new(),
        }
    }

    fn make_query_frame(seq: Vec<u8>, seg_masks: Vec<(usize, usize)>) -> QueryFrame {
        let mut aa_seq = Vec::with_capacity(seq.len() + 2);
        aa_seq.push(0);
        aa_seq.extend_from_slice(&seq);
        aa_seq.push(0);
        QueryFrame {
            frame: 1,
            aa_seq,
            aa_seq_nomask: None,
            aa_len: seq.len(),
            orig_len: seq.len(),
            seg_masks,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2803 (cited at
    // `x_context_karlin`)
    // The contexts of a query set with their Karlin blocks computed on a pool's threads
    // (`x_lookup_contexts_for_profile_in_pool`) are those of the one-thread computation: every
    // field, a context whose block fails (letters X only) included, for each matrix's bounds.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn contexts_in_pool_equal_the_one_thread_contexts() {
        let mut seed = 0x5eed_u64;
        let mut residue = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            1 + ((seed >> 33) % 22) as u8
        };
        let mut queries: Vec<Vec<QueryFrame>> = (0..37)
            .map(|i| {
                let seq: Vec<u8> = (0..(5 + i * 13)).map(|_| residue()).collect();
                vec![make_query_frame(seq, vec![(1, 3)])]
            })
            .collect();
        queries.push(vec![make_poly_a_query_frame()]);
        // Only negative scores: no Karlin block (NCBI's `is_valid = FALSE`).
        queries.push(vec![make_query_frame(vec![ncbistdaa::X; 20], vec![])]);
        for matrix in [
            ScoringMatrix::Blosum62,
            ScoringMatrix::Blosum45,
            ScoringMatrix::Pam30,
        ] {
            let one = x_lookup_contexts_for_profile(&queries, false, matrix, 3).unwrap();
            let pooled = crate::utils::threading::with_search_pool(4, "tblastn", |pool| {
                Ok(x_lookup_contexts_for_profile_in_pool(
                    &queries, false, matrix, 3, pool,
                ))
            })
            .unwrap()
            .unwrap();
            assert_eq!(one.len(), pooled.len());
            assert!(one.iter().any(|context| !context.is_valid));
            for (a, b) in one.iter().zip(&pooled) {
                assert_eq!(
                    (
                        a.q_idx,
                        a.f_idx,
                        a.frame,
                        a.aa_len,
                        a.orig_len,
                        a.frame_base
                    ),
                    (
                        b.q_idx,
                        b.f_idx,
                        b.frame,
                        b.aa_len,
                        b.orig_len,
                        b.frame_base
                    )
                );
                assert_eq!((&a.aa_seq, &a.aa_seq_nomask), (&b.aa_seq, &b.aa_seq_nomask));
                assert_eq!(a.is_valid, b.is_valid);
                assert_eq!(
                    format!("{:?}", a.karlin_params),
                    format!("{:?}", b.karlin_params)
                );
            }
        }
    }

    fn build_lookup_for_test(queries: Vec<Vec<QueryFrame>>, threshold: i32) -> BlastAaLookupTable {
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (lookup, _) = build_ncbi_lookup(&queries, threshold, &ungapped, true);
        lookup
    }

    fn prepare_lookup_for_test(queries: Vec<Vec<QueryFrame>>) -> PreparedLookupQuery {
        let ideal_params = blosum62_ideal_karlin_params();
        let std_comp = compute_std_aa_composition();
        prepare_lookup_query(
            &queries,
            ideal_params,
            &std_comp,
            LOOKUP_WORD_LENGTH,
            true,
            ScoringMatrix::Blosum62,
            (-4, 11),
        )
    }

    fn debruijn(alphabet_size: usize, order: usize) -> Vec<u8> {
        fn db(
            t: usize,
            p: usize,
            alphabet_size: usize,
            order: usize,
            a: &mut [usize],
            output: &mut Vec<u8>,
        ) {
            if t > order {
                if order % p == 0 {
                    for j in 1..=p {
                        output.push(a[j] as u8);
                    }
                }
                return;
            }

            a[t] = a[t - p];
            db(t + 1, p, alphabet_size, order, a, output);
            for residue in (a[t - p] + 1)..alphabet_size {
                a[t] = residue;
                db(t + 1, t, alphabet_size, order, a, output);
            }
        }

        let mut a = vec![0usize; alphabet_size * order + 1];
        let mut output = Vec::new();
        db(1, 1, alphabet_size, order, &mut a, &mut output);
        let wraparound = output
            .iter()
            .copied()
            .take(order.saturating_sub(1))
            .collect::<Vec<_>>();
        output.extend(wraparound);
        output
    }

    fn brute_force_lookup_count(subject_word: [u8; 3], threshold: i32) -> usize {
        let mut count = 0usize;
        for q0 in 0..BLASTAA_SIZE {
            for q1 in 0..BLASTAA_SIZE {
                for q2 in 0..BLASTAA_SIZE {
                    let score = blosum62_score(q0 as u8, subject_word[0])
                        + blosum62_score(q1 as u8, subject_word[1])
                        + blosum62_score(q2 as u8, subject_word[2]);
                    if score >= threshold
                        || (q0 as u8 == subject_word[0]
                            && q1 as u8 == subject_word[1]
                            && q2 as u8 == subject_word[2])
                    {
                        count += 1;
                    }
                }
            }
        }
        count
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-314
    // ```c
    // info->contexts[i].query_offset = new_offsets[i];
    // distance = new_offsets[i+1] - new_offsets[i];
    // info->contexts[i].query_length = distance ? distance-1 : 0;
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
    // ```c
    // offset += length + 1;
    // frame_offsets[context+1] = offset;
    // ```
    #[test]
    fn test_prepare_lookup_query_shares_single_boundary_sentinel_and_context_offsets() {
        let prepared = prepare_lookup_for_test(vec![
            vec![make_query_frame(
                vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                Vec::new(),
            )],
            vec![make_query_frame(
                vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
                Vec::new(),
            )],
        ]);

        assert_eq!(
            prepared.concat_query,
            vec![
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                0,
                ncbistdaa::C,
                ncbistdaa::C,
                ncbistdaa::C,
                0,
            ]
        );
        assert_eq!(prepared.frame_bases, vec![0, 4]);
        assert_eq!(prepared.lookup_locations, vec![(0, 2), (4, 6)]);
        assert_eq!(prepared.contexts[0].frame_base, 0);
        assert_eq!(prepared.contexts[1].frame_base, 4);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1045-1056
    // ```c
    // start_offset = query_info->contexts[context].query_offset;
    // end_offset = query_info->contexts[context].query_length + start_offset - 1;
    // if (mask_loc == NULL || mask_loc->seqloc_array[context] == NULL) {
    //     tail = BlastSeqLocNew(tail ? &tail : complement_mask,
    //                           start_offset, end_offset);
    //     continue;
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1067-1114
    // ```c
    // filter_start = start_offset + seq_range->left;
    // filter_end = start_offset + seq_range->right;
    // ...
    // tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
    // ```
    #[test]
    fn test_prepare_lookup_query_converts_seg_masks_to_absolute_inclusive_segments() {
        let prepared = prepare_lookup_for_test(vec![
            vec![make_query_frame(
                vec![
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                ],
                vec![(2, 4)],
            )],
            vec![make_query_frame(
                vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
                Vec::new(),
            )],
        ]);

        assert_eq!(prepared.lookup_locations, vec![(0, 1), (4, 5), (7, 9)]);
        assert_eq!(prepared.skipped_seg_mask, 4);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2789
    // ```c
    // loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    // if (loop_status) {
    //     contexts[context].is_valid = FALSE;
    //     sbp->sfp[context] = Blast_ScoreFreqFree(sbp->sfp[context]);
    //     sbp->kbp_std[context] = Blast_KarlinBlkFree(sbp->kbp_std[context]);
    //     continue;
    // }
    // ```
    #[test]
    fn test_prepare_lookup_query_marks_failed_ungapped_ka_context_invalid() {
        let prepared = prepare_lookup_for_test(vec![vec![make_query_frame(
            vec![ncbistdaa::X, ncbistdaa::X, ncbistdaa::X],
            Vec::new(),
        )]]);

        assert_eq!(prepared.contexts.len(), 1);
        assert!(!prepared.contexts[0].is_valid);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:246-250
    // ```c
    // return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-127
    // ```c
    // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
    // (*seq_blk)->length = length;
    // ```
    #[test]
    fn test_build_ncbi_lookup_query_length_matches_ncbi_sequence_length() {
        let lookup = build_lookup_for_test(
            vec![
                vec![make_query_frame(
                    vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                    Vec::new(),
                )],
                vec![make_query_frame(
                    vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
                    Vec::new(),
                )],
            ],
            0,
        );

        assert_eq!(lookup.query_length, 7);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:102-128
    // ```c
    // seq = query->sequence + from;
    // word_target = seq + lut_word_length;
    // ...
    // BlastLookupAddWordHit(backbone, lut_word_length, charsize,
    //                       seq - lut_word_length, offset - lut_word_length);
    // ```
    #[test]
    fn test_build_ncbi_lookup_multiple_queries_do_not_cross_shared_boundary_sentinel() {
        let lookup = build_lookup_for_test(
            vec![
                vec![make_query_frame(
                    vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                    Vec::new(),
                )],
                vec![make_query_frame(
                    vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                    Vec::new(),
                )],
            ],
            0,
        );
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.get_hits(aaa_index), &[0, 4]);
        assert_eq!(lookup.get_context_idx(0), 0);
        assert_eq!(lookup.get_context_idx(4), 1);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:340-357
    // ```c
    // if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
    //     dest = bbc[i].payload.entries;
    // else {
    //     bbc[i].payload.overflow_cursor = overflow_cursor;
    //     dest = (Int4 *) lookup->overflow + overflow_cursor;
    // }
    // for (j = 0; j < lookup->thin_backbone[i][1]; j++)
    //     dest[j] = lookup->thin_backbone[i][j + 2];
    // ```
    #[test]
    fn test_build_ncbi_lookup_preserves_overflow_query_offset_append_order() {
        let lookup = build_lookup_for_test(vec![vec![make_poly_a_query_frame()]], 0);
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.backbone[aaa_index].num_used, 4);
        assert_eq!(&lookup.overflow[..4], &[0, 1, 2, 3]);
        assert_eq!(lookup.get_hits(aaa_index), &[0, 1, 2, 3]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:504-543
    // ```c
    // if (lookup->threshold == 0 || score < lookup->threshold) {
    //     for (i = 0; i < offset_list[1]; i++) {
    //         BlastLookupAddWordHit(lookup->thin_backbone, lookup->word_length,
    //                               lookup->charsize, w,
    //                               query_bias + offset_list[i + 2]);
    //     }
    // }
    // ...
    // for (j = 0; j < offset_list[1]; j++) {
    //     BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
    //                           charsize, subject_word,
    //                           query_bias + offset_list[j + 2]);
    // }
    // ```
    #[test]
    fn test_build_ncbi_lookup_preserves_neighbor_query_offset_append_order() {
        let lookup = build_lookup_for_test(vec![vec![make_poly_a_query_frame()]], 11);
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.backbone[aaa_index].num_used, 4);
        assert_eq!(&lookup.overflow[..4], &[0, 1, 2, 3]);
        assert_eq!(lookup.get_hits(aaa_index), &[0, 1, 2, 3]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:96-128
    // ```c
    // if (word_length > to - from + 1)
    //     continue;
    // ...
    // if (*seq & invalid_mask)
    //     word_target = seq + lut_word_length + 1;
    // ```
    #[test]
    fn test_build_ncbi_lookup_skips_words_crossing_invalid_residue() {
        let lookup = build_lookup_for_test(
            vec![vec![make_query_frame(
                vec![
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    32,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                ],
                Vec::new(),
            )]],
            0,
        );
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.get_hits(aaa_index), &[0, 4]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:92-100
    // ```c
    // for (loc = locations; loc; loc = loc->next) {
    //     Int4 from = loc->ssr->left;
    //     Int4 to = loc->ssr->right;
    //     if (word_length > to - from + 1)
    //         continue;
    // ```
    #[test]
    fn test_build_ncbi_lookup_respects_interval_boundaries() {
        let lookup = build_lookup_for_test(
            vec![vec![make_query_frame(
                vec![
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                ],
                vec![(3, 4)],
            )]],
            0,
        );
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.get_hits(aaa_index), &[0]);
    }

    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:222-244 DebruijnSequenceTest
    // ```c
    // for(i=0;i<lookup->backbone_size;i++) {
    //   Int4 num_used = ((AaLookupSmallboneCell *)(lookup->thick_backbone))[i].num_used;
    //   if (((i & 0x1F) >= BLASTAA_SIZE) || (((i & 0x3E0) >> 5) >= BLASTAA_SIZE)
    //       || ((((i & 0x7C00) >> 10) >= BLASTAA_SIZE)))
    //       BOOST_REQUIRE_EQUAL(num_used, 0);
    //   else
    //       BOOST_REQUIRE_EQUAL(num_used, 1);
    // }
    // ```
    #[test]
    fn test_build_ncbi_lookup_debruijn_exact_words_match_ncbi_counts() {
        let seq = debruijn(BLASTAA_SIZE, LOOKUP_WORD_LENGTH);
        let lookup =
            build_lookup_for_test(vec![vec![make_query_frame(seq.clone(), Vec::new())]], 0);
        let charsize = lookup.charsize as usize;
        let residue_mask = (1usize << charsize) - 1;

        for idx in 0..lookup.backbone.len() {
            let w0 = (idx >> (2 * charsize)) & residue_mask;
            let w1 = (idx >> charsize) & residue_mask;
            let w2 = idx & residue_mask;
            if w0 < BLASTAA_SIZE && w1 < BLASTAA_SIZE && w2 < BLASTAA_SIZE {
                assert_eq!(lookup.backbone[idx].num_used, 1, "valid idx={idx}");
            } else {
                assert_eq!(lookup.backbone[idx].num_used, 0, "invalid idx={idx}");
            }
        }
    }

    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:253-277 NeighboringWordsTest
    // ```c
    // for each possible 3-mer:
    //   find its neighbors by brute force
    //   if (score >= 11) || exact-match
    //       count++;
    //   ensure that the number of neighbors matches the lookup table
    // ```
    #[test]
    fn test_build_ncbi_lookup_representative_neighbor_counts_match_bruteforce() {
        let seq = debruijn(BLASTAA_SIZE, LOOKUP_WORD_LENGTH);
        let lookup = build_lookup_for_test(vec![vec![make_query_frame(seq, Vec::new())]], 11);
        let representative_words = [
            [ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
            [ncbistdaa::W, ncbistdaa::W, ncbistdaa::W],
            [ncbistdaa::A, ncbistdaa::R, ncbistdaa::N],
        ];

        for subject_word in representative_words {
            let index = encode_kmer_3(
                subject_word[0] as usize,
                subject_word[1] as usize,
                subject_word[2] as usize,
                lookup.charsize as usize,
            ) & lookup.mask as usize;
            assert_eq!(
                lookup.backbone[index].num_used as usize,
                brute_force_lookup_count(subject_word, 11),
                "subject_word={subject_word:?}"
            );
        }
    }
}
