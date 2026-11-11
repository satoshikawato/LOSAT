//! One concatenated six-context BLASTX lookup and diagonal state.
use super::{
    args::ResolvedOptions, parameters::ContextParameters, query_setup::PreparedQueryBatch,
};
use crate::algorithm::blastp::extension::{
    extend_one_hit_blosum62, extend_one_hit_of, extend_two_hit_blosum62, extend_two_hit_table,
    BlastpOneHitUngappedResult, BlastpTwoHitUngappedResult,
};
use crate::algorithm::tblastx::blast_aascan::s_blast_aa_scan_subject_one_range;
use crate::algorithm::tblastx::blast_aascan::x_blast_aa_scan_subject_cells;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
    x_bx_build_compressed_lookup, x_bx_clut_pool, x_give_back_flat, x_give_back_wide_flat,
    XBxCompressedLookup, X_BLASTX,
};
use crate::algorithm::tblastx::lookup::compressed::x_compressed_scaled_threshold;
use crate::algorithm::tblastx::lookup::{
    build_ncbi_lookup_from_prepared, BackboneCell16, BlastAaLookupTable, QueryContext, XBxTable,
};
use crate::config::ScoringMatrix;
use anyhow::{ensure, Result};

/// No NCBI counterpart: the matrix of a search's seed and gapped stages as a type (item SC6 of
/// line NP), decided once before the search, so that BLOSUM62's instances ([`Blosum62Bx`]) call
/// today's functions (`extend_two_hit_blosum62`, `extend_one_hit_blosum62` and the stages written
/// for them) and keep their code, and another matrix ([`MatrixBx`]) gets instances of its own. NCBI
/// reads the matrix through `matrix->data` in the same loops (`aa_ungapped.c`
/// `s_BlastAaExtendTwoHit`, `s_BlastAaExtendOneHit`); the scores are the same.
pub(crate) trait BxScoring: Copy {
    /// Whether this is BLOSUM62's type: its stages are today's functions.
    const BLOSUM62: bool;
    /// The matrix.
    fn matrix(self) -> ScoringMatrix;
    /// The matrix's byte table of the two-hit walks (`x_ungapped_fast::matrix_table`), if any.
    fn table(self) -> Option<&'static crate::algorithm::blastp::x_ungapped_fast::MatrixTable>;
}

/// BLOSUM62's [`BxScoring`]: the default search's functions.
#[derive(Clone, Copy)]
pub(crate) struct Blosum62Bx;

impl BxScoring for Blosum62Bx {
    const BLOSUM62: bool = true;
    #[inline(always)]
    fn matrix(self) -> ScoringMatrix {
        ScoringMatrix::Blosum62
    }
    #[inline(always)]
    fn table(self) -> Option<&'static crate::algorithm::blastp::x_ungapped_fast::MatrixTable> {
        None
    }
}

/// Another matrix's [`BxScoring`]: BLASTP's extensions of its other matrices (item SC2b), the
/// two-hit walks on the matrix's byte table (`extend_two_hit_table`, the reference where it has
/// none) and the reference one-hit walk (`extend_one_hit_of`).
#[derive(Clone, Copy)]
pub(crate) struct MatrixBx {
    matrix: ScoringMatrix,
    table: Option<&'static crate::algorithm::blastp::x_ungapped_fast::MatrixTable>,
}

impl MatrixBx {
    pub(crate) fn of(matrix: ScoringMatrix) -> Self {
        debug_assert_ne!(matrix, ScoringMatrix::Blosum62);
        Self {
            matrix,
            table: crate::algorithm::blastp::x_ungapped_fast::matrix_table(matrix),
        }
    }
}

impl BxScoring for MatrixBx {
    const BLOSUM62: bool = false;
    #[inline(always)]
    fn matrix(self) -> ScoringMatrix {
        self.matrix
    }
    #[inline(always)]
    fn table(self) -> Option<&'static crate::algorithm::blastp::x_ungapped_fast::MatrixTable> {
        self.table
    }
}

// The extensions are free functions, not the trait's methods: rustc counts the inline methods of
// a trait impl as reachable from other crates (it cannot tell which impls generic code uses) and
// exports what they call, which would call BLOSUM62's `x_ungapped_fast::one_hit` through the GOT
// in today's loops.
/// `s_BlastAaExtendTwoHit` with the matrix of `scoring`: BLOSUM62's `extend_two_hit_blosum62`, or
/// the matrix's byte-table walks (`extend_two_hit_table`).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(super) fn bx_two_hit<S: BxScoring>(
    scoring: S,
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpTwoHitUngappedResult> {
    if S::BLOSUM62 {
        extend_two_hit_blosum62(
            query,
            subject,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    } else {
        extend_two_hit_table(
            (scoring.matrix(), scoring.table()),
            query,
            subject,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    }
}

/// `s_BlastAaExtendOneHit` with the matrix of `scoring`: BLOSUM62's `extend_one_hit_blosum62`, or
/// the reference walk of the matrix (`extend_one_hit_of`).
#[inline(always)]
pub(super) fn bx_one_hit<S: BxScoring>(
    scoring: S,
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<BlastpOneHitUngappedResult> {
    if S::BLOSUM62 {
        extend_one_hit_blosum62(query, subject, q_off, s_off, dropoff, word_size)
    } else {
        extend_one_hit_of(
            scoring.matrix(),
            query,
            subject,
            q_off,
            s_off,
            dropoff,
            word_size,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1040-1041
// ```c
//     aux_struct->offset_pairs =
//       (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
// ```
// The seed loop on a per-thread offset-pair buffer, read in place (a child module, so that it uses
// the diagonal table's fields).
#[path = "x_seed_scan.rs"]
pub mod x_seed_scan;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-511
// ```c
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
//         hits = scansub(lookup_wrap, subject,
//                                   offset_pairs, array_size, scan_range);
//
//         totalhits += hits;
//         /* for each hit, */
//         for (i = 0; i < hits; ++i) {
// ```
// EXPERIMENT (LOSAT_X_BXDIAG): the two-hit seed loop with the scan and the diagonal test fused
// (a child module, so that it uses the diagonal table's fields).
#[path = "x_seed_diag.rs"]
pub mod x_seed_diag;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155
// `s_BlastAaExtendTwoHit`.
// EXPERIMENT (LOSAT_X_BXUNGAPSENT): the two-hit extension of the seed loops with the walks that
// the sentinels stop.
#[path = "x_ungap_sent.rs"]
pub(crate) mod x_ungap_sent;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
// ```c++
//         Int4 diag_array_length;
//
//         diag_table= (BLAST_DiagTable*) calloc(1, sizeof(BLAST_DiagTable));
//
//         if (diag_table)
//         {
//                 diag_array_length = 1;
//                 /* What power of 2 is just longer than the query? */
//                 while (diag_array_length < (qlen+window_size))
//                 {
//                         diag_array_length = diag_array_length << 1;
//                 }
//                 /* These are used in the word finders to shift and mask
//                 rather than dividing and taking the remainder. */
//                 diag_table->diag_array_length = diag_array_length;
//                 diag_table->diag_mask = diag_array_length-1;
//                 diag_table->multiple_hits = multiple_hits;
//                 diag_table->offset = window_size;
//                 diag_table->window = window_size;
//         }
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
// ```c
// typedef struct DiagStruct {
//    signed int last_hit   : 31; /**< Offset of the last hit */
//    unsigned int flag      : 1 ; /**< Reset the next extension? */
// } DiagStruct;
// ```
// The entries are NCBI's bit-field (`DiagStruct`: `last_hit` truncated to 31 bits on write and
// sign-extended on read), so a `s + diag_offset` past 2^30 wraps as in NCBI. That happens at a
// `-window_size` near 2^29 from the second subject on (the offset grows by length + window per
// subject until INT4_MAX / 4); NCBI 2.17.0 then drops that subject's two-hit seeds, and so does
// LOSAT. The sums and differences on the entries wrap (two's complement i32), as NCBI's do.
pub struct Diagonals {
    entries: Vec<DiagStruct>,
    offset: i32,
    window: i32,
}
impl Diagonals {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
    // ```c++
    //         Int4 diag_array_length;
    //
    //         diag_table= (BLAST_DiagTable*) calloc(1, sizeof(BLAST_DiagTable));
    //
    //         if (diag_table)
    //         {
    //                 diag_array_length = 1;
    //                 /* What power of 2 is just longer than the query? */
    //                 while (diag_array_length < (qlen+window_size))
    //                 {
    //                         diag_array_length = diag_array_length << 1;
    //                 }
    //                 /* These are used in the word finders to shift and mask
    //                 rather than dividing and taking the remainder. */
    //                 diag_table->diag_array_length = diag_array_length;
    //                 diag_table->diag_mask = diag_array_length-1;
    //                 diag_table->multiple_hits = multiple_hits;
    //                 diag_table->offset = window_size;
    //                 diag_table->window = window_size;
    //         }
    // ```
    pub fn new(query_length: usize, window: i32) -> Result<Self> {
        Ok(Self {
            entries: vec![DiagStruct::default(); Self::array_length(query_length, window)?],
            offset: window,
            window,
        })
    }
    /// `diag_array_length` above: the power of two just at least `qlen + window_size`.
    fn array_length(query_length: usize, window: i32) -> Result<usize> {
        let n = query_length
            .checked_add(usize::try_from(window)?)
            .and_then(|n| n.max(1).checked_next_power_of_two())
            .ok_or_else(|| anyhow::anyhow!("BLASTX diagonal length exceeds Int4"))?;
        ensure!(
            n <= i32::MAX as usize,
            "BLASTX diagonal length exceeds Int4"
        );
        Ok(n)
    }
    // No NCBI counterpart: buffer reuse only (LOSAT_X_BXDIAGKEEP, blastx/x_bx_mem.rs). The entries,
    // offset and window `new(query_length, window)` gives, in this table's own allocation, which
    // `reserve_exact` grows where it is short (a reallocation, as `x_lut_direct::x_zeroed` grows
    // its vectors: a free of a mapped table would raise glibc's mmap and trim thresholds). A
    // thread that searches batches of two lengths keeps one table instead of freeing one and
    // making the other in turn.
    /// Makes this table the one `new(query_length, window)` makes.
    #[cfg(feature = "parallel")]
    pub(crate) fn x_remake(&mut self, query_length: usize, window: i32) -> Result<()> {
        let n = Self::array_length(query_length, window)?;
        self.entries.clear();
        self.entries.reserve_exact(n);
        self.entries.resize(n, DiagStruct::default());
        self.offset = window;
        self.window = window;
        Ok(())
    }
    // No NCBI counterpart: buffer reuse only (LOSAT_X_BXDIAGKEEP, blastx/x_bx_mem.rs). The same
    // window and the same power-of-two length give the same mask, so the table's diagonals collide
    // as those of a table `new(query_length, window)` makes.
    /// Whether this table has the length and window `new(query_length, window)` gives.
    #[cfg(feature = "parallel")]
    pub(crate) fn x_fits(&self, query_length: usize, window: i32) -> bool {
        self.window == window
            && Self::array_length(query_length, window).ok() == Some(self.entries.len())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-175
    // ```c++
    // Blast_ExtendWordExit(Blast_ExtendWord * ewp, Int4 subject_length)
    // {
    //     if (!ewp)
    //         return -1;
    //
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    //         if (ewp->hash_table->offset >= INT4_MAX / 4) {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
    // ```c++
    // static Int4 s_BlastDiagClear(BLAST_DiagTable * diag)
    // {
    //     Int4 i, n;
    //     DiagStruct *diag_struct_array;
    //
    //     if (diag == NULL)
    //         return 0;
    //
    //     n = diag->diag_array_length;
    //
    //     diag->offset = diag->window;
    //
    //     diag_struct_array = diag->hit_level_array;
    //
    //     for (i = 0; i < n; i++) {
    //         diag_struct_array[i].flag = 0;
    //         diag_struct_array[i].last_hit = -diag->window;
    //         if (diag->hit_len_array) diag->hit_len_array[i] = 0;
    //     }
    //     return 0;
    // }
    // ```
    pub fn finish_subject(&mut self, length: usize) -> Result<()> {
        if self.offset >= i32::MAX / 4 {
            self.offset = self.window;
            self.entries.fill(DiagStruct::clear(self.window));
        } else {
            self.offset = self
                .offset
                .checked_add(i32::try_from(length)?)
                .and_then(|v| v.checked_add(self.window))
                .ok_or_else(|| anyhow::anyhow!("BLASTX diagonal offset exceeds Int4"))?;
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:273-313
// ```c++
//  */
// static int score_compare_match(const void *v1, const void *v2)
// {
//     BlastInitHSP *h1, *h2;
//     int result = 0;
//
//     h1 = (BlastInitHSP *) v1;
//     h2 = (BlastInitHSP *) v2;
//
//     /* Check if ungapped_data substructures are initialized. If not, move
//        those array elements to the end. In reality this should never happen. */
//     if (h1->ungapped_data == NULL && h2->ungapped_data == NULL)
//         return 0;
//     else if (h1->ungapped_data == NULL)
//         return 1;
//     else if (h2->ungapped_data == NULL)
//         return -1;
//
//     if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
//                                  h1->ungapped_data->score)) &&
//         0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
//                                  h2->ungapped_data->s_start)) &&
//         0 == (result = BLAST_CMP(h2->ungapped_data->length,
//                                  h1->ungapped_data->length)) &&
//         0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
//                                  h2->ungapped_data->q_start))) {
//         result = BLAST_CMP(h2->ungapped_data->length,
//                            h1->ungapped_data->length);
//     }
//
//     return result;
// }
//
// void Blast_InitHitListSortByScore(BlastInitHitList * init_hitlist)
// {
//     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
//           sizeof(BlastInitHSP), score_compare_match);
// }
//
// Boolean Blast_InitHitListIsSortedByScore(BlastInitHitList * init_hitlist)
// {
// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InitHsp {
    pub q_seed: i32,
    pub s_seed: i32,
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:484-496
// ```c++
//         scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
//         wordsize = lookup->word_length;
//         use_pssm = lookup->use_pssm;
//     }
//     else {
//         BlastCompressedAaLookupTable *lookup =
//                         (BlastCompressedAaLookupTable *)(lookup_wrap->lut);
//         scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
//         wordsize = lookup->word_length;
//     }
//
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
// ```
// EXPERIMENT (LOSAT_X_BXSCANBYTE): the standard table carries the presence bytes of the direct
// loop's scan, built the first time that loop scans with it (x_seed_scan.rs).
#[expect(
    clippy::large_enum_variant,
    reason = "boxing a variant changes the layout and the memory; that needs an RSS and time measurement"
)]
pub(crate) enum Lookup {
    Standard(BlastAaLookupTable, x_seed_scan::XKeptPresence),
    // EXPERIMENT (LOSAT_X_BPCLUT): the compressed table as the reference builds it or flat
    // (`x_clut_direct.rs`), one variant so that the standard table's loops keep their match.
    Compressed(XBxCompressedLookup),
}

// No NCBI counterpart: buffer reuse only (the pooled build, tblastx/lookup/x_lut_direct.rs): a
// dropped table's vectors are kept for the next build instead of freed, where the table was built
// pooled (not in diagnostics mode, whose reference builder makes the tables, nor with
// LOSAT_X_BXPOOLN1 in a search on one thread) and no panic unwinds.
impl Drop for Lookup {
    fn drop(&mut self) {
        match self {
            Self::Standard(table, kept) => {
                if kept.pooled && !std::thread::panicking() {
                    crate::algorithm::tblastx::lookup::x_give_back_pooled(table);
                }
            }
            Self::Compressed(table) => x_give_back_compressed(table),
        }
    }
}

// No NCBI counterpart: the BLASTX lookup table of a query batch, either width of the standard
// table's hits. `Wide` is `Lookup` as before 16-bit hits (the table of `Int4` hits as NCBI's, or the
// compressed table), which the direct loop and its scans read. `Narrow` is the standard table with
// its hits stored in 16 bits (tblastx/lookup/backbone.rs `BackboneCell16`), built where the fused
// loop scans the table (`x_seed_diag::x_narrow_cutoff`: the query batch fits 16 bits and the table
// has at least LOSAT_X_BXDIAGMIN hits per used cell). It is a type of its own, so that no loop of
// `Int4` hits can be given one: the fused loop is generic over the table's cell and scans either,
// and the reference scan reads either.
#[expect(
    clippy::large_enum_variant,
    reason = "boxing a variant changes the layout and the memory; that needs an RSS and time measurement"
)]
pub(crate) enum BxLookup {
    Wide(Lookup),
    Narrow(NarrowLookup),
}

/// The standard lookup table with 16-bit hits and its companion (`XKeptPresence`: the build
/// number, the context table, whether the table goes back to the pool); see `BxLookup`.
pub struct NarrowLookup {
    pub(crate) table: BlastAaLookupTable<BackboneCell16>,
    pub(crate) kept: x_seed_scan::XKeptPresence,
}

impl BxLookup {
    /// [`Lookup::x_new_as`] with the cutoff of the search's options (`x_seed_diag::x_narrow_cutoff`).
    pub(crate) fn x_new_in_parts(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
        parts: usize,
    ) -> Result<Self> {
        // The search's matrix chooses the builds once per query batch: BLOSUM62's are
        // `x_new_as`'s, another matrix's `x_new_of`'s (the same tables with its scores).
        if options.scoring != ScoringMatrix::Blosum62 {
            return Lookup::x_new_of(
                batch,
                parameters,
                options,
                parts,
                x_seed_diag::x_narrow_cutoff(options),
            );
        }
        Lookup::x_new_as(
            batch,
            parameters,
            options,
            parts,
            x_seed_diag::x_narrow_cutoff(options),
        )
    }
}

// No NCBI counterpart: buffer reuse only (as the drop of `Lookup`): the vectors of a table of
// 16-bit hits go back to the pool when the table was built pooled and no panic unwinds.
impl Drop for NarrowLookup {
    fn drop(&mut self) {
        if self.kept.pooled && !std::thread::panicking() {
            crate::algorithm::tblastx::lookup::x_give_back_pooled16(&mut self.table);
        }
    }
}

impl NarrowLookup {
    /// [`Lookup::scan`] for the table of 16-bit hits (the reference scan); `ranges` are the scan
    /// ranges of a `-subject` with lower-case letters (`x_soft_subject.rs`), or the whole subject
    /// (`None`).
    pub(super) fn x_scan(
        &self,
        subject: &[u8],
        length: usize,
        ranges: Option<&[(i32, i32)]>,
        window: i32,
        mut visit: impl FnMut(&[(i32, i32)]),
    ) {
        let word = self.table.word_length;
        let whole = [(0, length as i32)];
        let ranges = ranges.unwrap_or(&whole);
        let mut range = [0, ranges[0].0, ranges[0].1 - word];
        if window > 0 && range[2] < range[1] {
            range[2] = range[1];
        }
        let mut pad = Vec::new();
        let subject = x_scan_bytes(subject, range[2], word, &mut pad);
        let mut pairs = vec![OffsetPair::default(); 4096 + self.table.longest_chain as usize];
        let capacity = pairs.len() as i32;
        let mut seeds = Vec::new();
        while range[1] <= range[2] {
            seeds.clear();
            let n = x_blast_aa_scan_subject_cells(
                &self.table,
                subject,
                ranges,
                &mut pairs,
                capacity,
                &mut range,
            );
            for pair in &pairs[..n as usize] {
                seeds.push((pair.q_off as i32, pair.s_off as i32));
            }
            visit(&seeds);
        }
    }
}

// No NCBI counterpart: buffer reuse only (LOSAT_X_BXCLUTPOOL, tblastx/lookup/x_clut_direct.rs): a
// dropped flat table's vectors are kept for the next build. Not inlined, so that a standard table's
// drop keeps its code.
#[cold]
#[inline(never)]
fn x_give_back_compressed(table: &mut XBxCompressedLookup) {
    if !x_bx_clut_pool() || std::thread::panicking() {
        return;
    }
    match table {
        XBxCompressedLookup::Reference(_) => {}
        XBxCompressedLookup::Flat(flat) => x_give_back_flat(flat),
        XBxCompressedLookup::Wide(flat) => x_give_back_wide_flat(flat),
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// ```c
//     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                     query, locations);
//     s_CompressedLookupFinalize(lookup);
// ```
/// EXPERIMENT (LOSAT_X_BPCLUT): the compressed table of a query batch (`-word_size 5`), flat as
/// TBLASTN's with `u16` query offsets (`i32` past `X_U16_LETTERS` letters;
/// `x_bx_build_compressed_lookup`, on vectors kept between builds with `LOSAT_X_BXCLUTPOOL`), on one
/// thread. One table for the batch, so its one
/// range is all the segments. Words that fill NCBI's overflow banks are refused, as BLASTP's
/// (NCBI writes past its banks). Cold and not inlined, so that the standard table's build keeps
/// its code.
#[cold]
#[inline(never)]
fn x_compressed_lookup(query: &[u8], segments: &[(i32, i32)], threshold: f64) -> Result<Lookup> {
    x_bx_build_compressed_lookup(
        5,
        x_compressed_scaled_threshold(threshold),
        query,
        segments,
        &[0..segments.len()],
    )
    .ok_or_else(|| anyhow::anyhow!("unsupported BLASTX compressed lookup"))?
    .map(Lookup::Compressed)
    .map_err(|()| {
        anyhow::anyhow!(
            "-threshold {} with the compressed lookup table of -word_size 5 is not supported by LOSAT's BLASTX: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap",
            threshold
        )
    })
}

/// No NCBI counterpart: the word length of either form of the compressed table. Not inlined, so
/// that the standard table's loops keep their code.
#[inline(never)]
pub(super) fn x_compressed_word_length(table: &XBxCompressedLookup) -> i32 {
    match table {
        XBxCompressedLookup::Reference(table) => table.word_length,
        XBxCompressedLookup::Flat(table) => table.word_length,
        XBxCompressedLookup::Wide(table) => table.word_length,
    }
}

/// No NCBI counterpart: `scan_subject` of either form of the compressed table, which give the same
/// pairs and `s_range`. Not inlined, so that the standard table's loops keep their code.
#[inline(never)]
pub(super) fn x_compressed_scan(
    table: &XBxCompressedLookup,
    subject: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    match table {
        XBxCompressedLookup::Reference(table) => {
            table.scan_subject(subject, seq_ranges, offset_pairs, array_size, s_range)
        }
        XBxCompressedLookup::Flat(table) => {
            table.scan_subject::<X_BLASTX>(subject, seq_ranges, offset_pairs, array_size, s_range)
        }
        XBxCompressedLookup::Wide(table) => {
            table.scan_subject::<X_BLASTX>(subject, seq_ranges, offset_pairs, array_size, s_range)
        }
    }
}

/// The contexts of the lookup query of `batch` (`x_new_as`, `x_new_of`): each context's letters
/// between sentinels, its frame and its ungapped block.
#[inline(always)]
fn lookup_contexts(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    query: &[u8],
) -> Vec<QueryContext> {
    batch
        .contexts
        .iter()
        .zip(parameters)
        .enumerate()
        .map(|(i, (c, p))| {
            let mut aa = vec![0];
            aa.extend_from_slice(&query[c.offset..c.offset + c.length]);
            aa.push(0);
            QueryContext {
                q_idx: c.query_index as u32,
                f_idx: (i % 6) as u8,
                frame: c.frame,
                aa_seq: aa,
                aa_seq_nomask: None,
                aa_len: c.length,
                orig_len: batch.original_lengths[c.query_index],
                frame_base: c.offset as i32,
                is_valid: p.valid,
                karlin_params: p.ungapped,
            }
        })
        .collect()
}

impl Lookup {
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1349-1351
    // ```c++
    //
    //     /* index the query and finish up */
    //
    // ```
    // Only the tests build a table this way; the searches use `x_new_as`.
    #[cfg(test)]
    pub(crate) fn new(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
    ) -> Result<Self> {
        match Self::x_new_as(batch, parameters, options, 1, None)? {
            BxLookup::Wide(lookup) => Ok(lookup),
            BxLookup::Narrow(_) => unreachable!("a build without a cutoff takes Int4 hits"),
        }
    }

    /// EXPERIMENT (LOSAT_X_BXBUILDPAR): [`Lookup::new`] with a pooled standard table's passes on
    /// `parts` tasks of the pool it runs on (`tblastx/lookup/x_lut_direct.rs`
    /// `x_partitioned_passes`); the table is the same. The standard table takes 16-bit hits (a
    /// [`NarrowLookup`]) where the lookup query fits and its hits per used cell, in hundredths,
    /// reach `narrow_cutoff` (`None`: never).
    pub(crate) fn x_new_as(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
        parts: usize,
        narrow_cutoff: Option<u64>,
    ) -> Result<BxLookup> {
        let query = batch.sequence_start[1..].to_vec();
        if options.word_size == 5 {
            return x_compressed_lookup(&query, &batch.lookup_segments, options.threshold)
                .map(BxLookup::Wide);
        }
        ensure!(
            options.word_size == 3,
            "unsupported BLASTX lookup word size"
        );
        let contexts = lookup_contexts(batch, parameters, &query);
        // Dispatch point (LOSAT_X_BXPOOLN1, in the strict set): a search on one thread builds its
        // tables without the pool (tblastx/lookup/x_lut_direct.rs). No NCBI counterpart.
        let pooled = crate::algorithm::tblastx::lookup::x_bx_pooled_for(
            usize::try_from(options.num_threads).unwrap_or(1),
        );
        let kept = x_seed_scan::XKeptPresence::x_built(pooled);
        Ok(
            match build_ncbi_lookup_from_prepared(
                query,
                batch.lookup_segments.clone(),
                contexts,
                // `lookup->threshold = (Int4)opt->threshold` (blast_aalookup.c:241), converted as
                // the x86-64 build does (`INT_MIN` past `Int4`).
                crate::core::blast_util::ncbi_int4_from_double(options.threshold),
                pooled,
                parts,
                narrow_cutoff,
            ) {
                XBxTable::Wide(table) => BxLookup::Wide(Self::Standard(table, kept)),
                XBxTable::Narrow(table) => BxLookup::Narrow(NarrowLookup { table, kept }),
            },
        )
    }
    /// [`Lookup::x_new_as`] with another matrix than BLOSUM62 (items SC6 and SC6b of line NP): the
    /// compressed table of `-word_size 5` with the matrix's scores (`x_tbn_build_compressed_lookup_of`,
    /// as TBLASTN's, `Int4` hits, without the pool of kept vectors), the standard one as BLOSUM62's
    /// (`build_ncbi_lookup_from_prepared_of`: the pool, `parts`, 16-bit hits from `narrow_cutoff`)
    /// with its scores. Not inlined, so that BLOSUM62's builds keep their code.
    #[cold]
    #[inline(never)]
    pub(crate) fn x_new_of(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
        parts: usize,
        narrow_cutoff: Option<u64>,
    ) -> Result<BxLookup> {
        use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
            x_tbn_build_compressed_lookup_of, XCompressedLookup,
        };
        let matrix = options.scoring;
        let query = batch.sequence_start[1..].to_vec();
        if options.word_size == 5 {
            let segments = &batch.lookup_segments;
            return x_tbn_build_compressed_lookup_of(
                matrix,
                5,
                x_compressed_scaled_threshold(options.threshold),
                &query,
                segments,
                &[0..segments.len()],
            )
            .ok_or_else(|| anyhow::anyhow!("unsupported BLASTX compressed lookup"))?
            .map(|table| {
                BxLookup::Wide(Lookup::Compressed(match table {
                    XCompressedLookup::Reference(table) => XBxCompressedLookup::Reference(table),
                    XCompressedLookup::Flat(table) => XBxCompressedLookup::Wide(table),
                }))
            })
            .map_err(|()| {
                anyhow::anyhow!(
                    "-threshold {} with the compressed lookup table of -word_size 5 is not supported by LOSAT's BLASTX: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap",
                    options.threshold
                )
            });
        }
        ensure!(
            options.word_size == 3,
            "unsupported BLASTX lookup word size"
        );
        let contexts = lookup_contexts(batch, parameters, &query);
        // The pool as `x_new_as` takes it (LOSAT_X_BXPOOLN1).
        let pooled = crate::algorithm::tblastx::lookup::x_bx_pooled_for(
            usize::try_from(options.num_threads).unwrap_or(1),
        );
        let kept = x_seed_scan::XKeptPresence::x_built(pooled);
        Ok(
            match crate::algorithm::tblastx::lookup::build_ncbi_lookup_from_prepared_of(
                matrix,
                query,
                batch.lookup_segments.clone(),
                contexts,
                // `lookup->threshold = (Int4)opt->threshold` (blast_aalookup.c:241), as `x_new_as`.
                crate::core::blast_util::ncbi_int4_from_double(options.threshold),
                pooled,
                parts,
                narrow_cutoff,
            ) {
                XBxTable::Wide(table) => BxLookup::Wide(Self::Standard(table, kept)),
                XBxTable::Narrow(table) => BxLookup::Narrow(NarrowLookup { table, kept }),
            },
        )
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
    // ```c++
    // Int4 GetOffsetArraySize(LookupTableWrap* lookup)
    // {
    //    Int4 offset_array_size;
    //
    //    switch (lookup->lut_type) {
    //    case eMBLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastMBLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastAaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eCompressedAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastCompressedAaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eSmallNaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastSmallNaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eNaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastNaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eNaHashLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastNaHashLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    default:
    //       offset_array_size = OFFSET_ARRAY_SIZE;
    //       break;
    //    }
    //    return offset_array_size;
    // ```
    pub fn capacity(&self) -> usize {
        4096 + match self {
            Self::Standard(l, _) => l.longest_chain as usize,
            Self::Compressed(XBxCompressedLookup::Reference(l)) => l.longest_chain as usize,
            Self::Compressed(XBxCompressedLookup::Flat(l)) => l.longest_chain as usize,
            Self::Compressed(XBxCompressedLookup::Wide(l)) => l.longest_chain as usize,
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
    // ```c++
    //         wordsize = lookup->word_length;
    //     }
    //
    //     scan_range[0] = 0;
    //     scan_range[1] = subject->seq_ranges[0].left;
    //     scan_range[2] = subject->seq_ranges[0].right - wordsize;
    //
    //     if (scan_range[2] < scan_range[1])
    //         scan_range[2] = scan_range[1];
    //
    //     while (scan_range[1] <= scan_range[2]) {
    //         /* scan the subject sequence for hits */
    //         hits = scansub(lookup_wrap, subject,
    //                                   offset_pairs, array_size, scan_range);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-769
    // ```c++
    //     scan_range[0] = 0;
    //     scan_range[1] = subject->seq_ranges[0].left;
    //     scan_range[2] = subject->seq_ranges[0].right - wordsize;
    //
    //     while (scan_range[1] <= scan_range[2]) {
    //         /* scan the subject sequence for hits */
    //         hits = scansub(lookup_wrap, subject,
    //                        offset_pairs, array_size, scan_range);
    // ```
    pub(crate) fn scan(
        &self,
        subject: &[u8],
        length: usize,
        window: i32,
        mut visit: impl FnMut(&[(i32, i32)]),
    ) {
        let word = match self {
            Self::Standard(l, _) => l.word_length,
            Self::Compressed(l) => x_compressed_word_length(l),
        };
        let ranges = [(0, length as i32)];
        let mut range = [0, 0, length as i32 - word];
        if window > 0 && range[2] < range[1] {
            range[2] = range[1];
        }
        let mut pad = Vec::new();
        let subject = x_scan_bytes(subject, range[2], word, &mut pad);
        let mut pairs = vec![OffsetPair::default(); self.capacity()];
        let capacity = pairs.len() as i32;
        let mut seeds = Vec::new();
        while range[1] <= range[2] {
            seeds.clear();
            let n = match self {
                Self::Standard(l, _) => {
                    s_blast_aa_scan_subject_one_range(l, subject, &mut pairs, capacity, &mut range)
                }
                Self::Compressed(l) => {
                    x_compressed_scan(l, subject, &ranges, &mut pairs, capacity, &mut range)
                }
            };
            for pair in &pairs[..n as usize] {
                seeds.push((pair.q_off as i32, pair.s_off as i32));
            }
            visit(&seeds);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:505-506
// ```c
//     if (scan_range[2] < scan_range[1])
//         scan_range[2] = scan_range[1];
// ```
// With a two-hit window the last scan position is at least 0, so the one position of a subject
// shorter than a word less one letter has a word that runs past the sentinel, and NCBI reads the
// bytes after it in its buffer (the next sequence of a database volume). No NCBI counterpart for
// the copy: LOSAT scans such a subject from a copy with sentinels after it, as the fused loop
// reads it (`blastp::x_seed_fused::x_scan_cells_bf`). One scan position extends no hit (the
// two-hit test needs two hits on a diagonal), and the entry such a hit leaves lies a window and
// the subject's length behind the next subject's offsets, so no hit and no report depends on
// those bytes.
/// The bytes a scan of positions up to `s_last` (each with its `word` letters) reads: `subject`
/// (the residues and the sentinel), or, when the last word runs past it, a copy in `pad` with
/// sentinels (0) after it.
pub(crate) fn x_scan_bytes<'a>(
    subject: &'a [u8],
    s_last: i32,
    word: i32,
    pad: &'a mut Vec<u8>,
) -> &'a [u8] {
    let read = (s_last + word).max(0) as usize;
    if read <= subject.len() {
        return subject;
    }
    pad.clear();
    pad.extend_from_slice(subject);
    pad.resize(read, 0);
    pad
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-607
// ```c++
//             diag_coord = (query_offset - subject_offset) & diag_mask;
//
//             /* If the reset bit is set, an extension just happened. */
//             if (diag_array[diag_coord].flag) {
//                 /* If we've already extended past this hit, skip it. */
//                 if ((Int4) (subject_offset + diag_offset) <
//                     diag_array[diag_coord].last_hit) {
//                     continue;
//                 }
//                 /* Otherwise, start a new hit. */
//                 else {
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     diag_array[diag_coord].flag = 0;
//                 }
//             }
//             /* If the reset bit is cleared, try to start an extension. */
//             else {
//                 /* find the distance to the last hit on this diagonal */
//                 last_hit = diag_array[diag_coord].last_hit - diag_offset;
//                 diff = subject_offset - last_hit;
//
//                 if (diff >= window) {
//                     /* We are beyond the window for this diagonal; start a
//                        new hit */
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     continue;
//                 }
//
//                 /* If the difference is less than the wordsize (i.e. last
//                    hit and this hit overlap), give up */
//
//                 if (diff < wordsize) {
//                     continue;
//                 }
//
//                 /* Extend this pair of hits. The extension to the left must
//                    reach the end of the first word in order for extension to
//                    the right to proceed.
//
//                    To use the cutoff and X-drop values appropriate for this
//                    extension, the query context must first be found */
//
//                 curr_context = BSearchContextInfo(query_offset, query_info);
//
//                 /* Check if the last hit hits current query. Because last_hit
//                    is never reset, it may contain a hit to the previous
//                    concatenated query */
//
//                 if (query_offset - diff <
//                     query_info->contexts[curr_context].query_offset) {
//
//                     /* there was no last hit for this diagnol; start a new hit */
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     continue;
//                 }
//
//                 cutoffs = word_params->cutoffs + curr_context;
//                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
//                                               last_hit + wordsize,
//                                               subject_offset, query_offset,
//                                               cutoffs->x_dropoff,
//                                               &hsp_q, &hsp_s,
//                                               &hsp_len, use_pssm,
//                                               wordsize, &right_extend,
//                                               &s_last_off);
//
//                 ++hits_extended;
//
//                 /* if the hsp meets the score threshold, report it */
//                 if (score >= cutoffs->cutoff_score)
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//                                      query_offset, subject_offset, hsp_len,
//                                      score);
//
//                 /* If an extension to the right happened, reset the last hit
//                    so that future hits to this diagonal must start over. */
//
//                 if (right_extend) {
//                     diag_array[diag_coord].flag = 1;
//                     diag_array[diag_coord].last_hit =
//                         s_last_off - (wordsize - 1) + diag_offset;
//                 }
//                 /* Otherwise, make the present hit into the previous hit for
//                    this diagonal */
//                 else {
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                 }
//             }                   /* end else */
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:774-800
// ```c++
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//             diag_coord = (subject_offset - query_offset) & diag_mask;
//             diff = subject_offset -
//                 (diag_array[diag_coord].last_hit - diag_offset);
//
//             /* do an extension, but only if we have not already extended this
//                far */
//             if (diff >= 0) {
//                 Int4 curr_context = BSearchContextInfo(query_offset,
//                                                        query_info);
//                 BlastUngappedCutoffs *cutoffs = word_params->cutoffs +
//                                                         curr_context;
//                 score = s_BlastAaExtendOneHit(matrix, subject, query,
//                                               subject_offset, query_offset,
//                                               cutoffs->x_dropoff,
//                                               &hsp_q, &hsp_s, &hsp_len,
//                                               wordsize, use_pssm, &s_last_off);
//
//                 /* if the hsp meets the score threshold, report it */
//                 if (score >= cutoffs->cutoff_score) {
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//                                      query_offset, subject_offset, hsp_len,
//                                      score);
//                 }
//                 diag_array[diag_coord].last_hit =
//                     s_last_off - (wordsize - 1) + diag_offset;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-515
// ```c++
//         scan_range[2] = scan_range[1];
//
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
//         hits = scansub(lookup_wrap, subject,
//                                   offset_pairs, array_size, scan_range);
//
//         totalhits += hits;
//         /* for each hit, */
//         for (i = 0; i < hits; ++i) {
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//
//             /* calculate the diagonal associated with this query-subject pair
//              */
//
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ...
//                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
// ...
//                 if (score >= cutoffs->cutoff_score)
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
// ...
//                 if (right_extend) {
//                     diag_array[diag_coord].flag = 1;
//                     diag_array[diag_coord].last_hit =
//                         s_last_off - (wordsize - 1) + diag_offset;
// ```
// word_finder is the original port of this loop for one subject; the diagnostic stages run it
// (they trace the seeds), and the unit tests of the faster loops compare with it.
/// The reference word finder (the loop of NCBI's BlastAaWordFinder_TwoHit /
/// _OneHit for one subject), kept verbatim.
pub(crate) fn word_finder(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    word_finder_over(
        Blosum62Bx,
        batch,
        parameters,
        options,
        subject,
        |each| {
            lookup.scan(
                subject,
                subject.len().saturating_sub(1),
                options.window_size,
                each,
            )
        },
        diagonals,
        trace,
    )
}

/// [`word_finder_bx`] with another matrix than BLOSUM62 (`MatrixBx`).
#[cold]
#[inline(never)]
pub(super) fn word_finder_bx_of(
    scoring: MatrixBx,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    match lookup {
        BxLookup::Wide(lookup) => word_finder_any(
            scoring, batch, parameters, options, subject, lookup, diagonals, trace,
        ),
        BxLookup::Narrow(lookup) => word_finder_narrow(
            scoring, batch, parameters, options, subject, lookup, diagonals, trace,
        ),
    }
}

/// [`word_finder`] with the matrix of `scoring`: BLOSUM62's is `word_finder` itself.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(super) fn word_finder_any<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    if S::BLOSUM62 {
        return word_finder(
            batch, parameters, options, subject, lookup, diagonals, trace,
        );
    }
    word_finder_over(
        scoring,
        batch,
        parameters,
        options,
        subject,
        |each| {
            lookup.scan(
                subject,
                subject.len().saturating_sub(1),
                options.window_size,
                each,
            )
        },
        diagonals,
        trace,
    )
}

/// [`word_finder`] for a lookup table of either width ([`BxLookup`]).
pub(super) fn word_finder_bx(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    match lookup {
        BxLookup::Wide(lookup) => word_finder(
            batch, parameters, options, subject, lookup, diagonals, trace,
        ),
        BxLookup::Narrow(lookup) => word_finder_narrow(
            Blosum62Bx, batch, parameters, options, subject, lookup, diagonals, trace,
        ),
    }
}

/// [`word_finder`] for the table of 16-bit hits, with the matrix of `scoring`.
#[allow(clippy::too_many_arguments)]
pub(super) fn word_finder_narrow<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &NarrowLookup,
    diagonals: &mut Diagonals,
    trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    word_finder_over(
        scoring,
        batch,
        parameters,
        options,
        subject,
        |each| {
            lookup.x_scan(
                subject,
                subject.len().saturating_sub(1),
                None,
                options.window_size,
                each,
            )
        },
        diagonals,
        trace,
    )
}

/// [`word_finder`] with the scan of the subject given as `scan`, which calls its argument with each
/// batch of the scanner's pairs: the whole subject (`Lookup::scan`), or the scan ranges of a
/// `-subject` with lower-case letters (`-lcase_masking`, `x_soft_subject.rs`).
pub(super) fn word_finder_over<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    scan: impl FnOnce(&mut dyn FnMut(&[(i32, i32)])),
    diagonals: &mut Diagonals,
    mut trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    let query = &batch.sequence_start[1..];
    let word = options.word_size;
    let mask = diagonals.entries.len() as i32 - 1;
    let mut hits = Vec::new();
    let mut each = |seeds: &[(i32, i32)]| {
        trace(seeds);
        for &(q, s) in seeds {
            // The context of a seed is only needed once the diagonal tests let it through (the
            // lookup has no side effect).
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560,783-784
            // ```c
            //                 curr_context = BSearchContextInfo(query_offset, query_info);
            // ...
            //                 Int4 curr_context = BSearchContextInfo(query_offset,
            //                                                        query_info);
            // ```
            // The two-hit path looks up the context (BSearchContextInfo) only after the
            // diagonal tests have let the pair through, and the one-hit path only when diff >=
            // 0. The same here.
            let context_of = |q: i32| {
                batch
                    .contexts
                    .partition_point(|c| c.offset as i32 <= q)
                    .saturating_sub(1)
            };
            let context;
            let index = if options.window_size == 0 {
                (s - q) & mask
            } else {
                (q - s) & mask
            } as usize;
            let entry = &mut diagonals.entries[index];
            let (u, end, extended) = if options.window_size == 0 {
                if s.wrapping_sub(entry.last_hit().wrapping_sub(diagonals.offset)) < 0 {
                    continue;
                }
                context = context_of(q);
                let p = &parameters[context];
                let Some(result) = bx_one_hit(
                    scoring,
                    query,
                    subject,
                    q as usize,
                    s as usize,
                    p.word_xdrop,
                    word as usize,
                ) else {
                    continue;
                };
                (result.ungapped_data, result.s_last_off, true)
            } else {
                if entry.flag() != 0 {
                    if s.wrapping_add(diagonals.offset) < entry.last_hit() {
                        continue;
                    }
                    entry.set_last_hit(s.wrapping_add(diagonals.offset));
                    entry.set_flag(0);
                    continue;
                }
                let previous = entry.last_hit().wrapping_sub(diagonals.offset);
                let diff = s.wrapping_sub(previous);
                if diff >= options.window_size {
                    entry.set_last_hit(s.wrapping_add(diagonals.offset));
                    continue;
                }
                if diff < word {
                    continue;
                }
                context = context_of(q);
                let p = &parameters[context];
                if q - diff < batch.contexts[context].offset as i32 {
                    entry.set_last_hit(s.wrapping_add(diagonals.offset));
                    continue;
                }
                let Some(result) = bx_two_hit(
                    scoring,
                    query,
                    subject,
                    (previous + word) as usize,
                    s as usize,
                    q as usize,
                    p.word_xdrop,
                    word as usize,
                ) else {
                    continue;
                };
                (result.ungapped_data, result.s_last_off, result.right_extend)
            };
            let p = &parameters[context];
            if u.score >= p.word_cutoff {
                hits.push(InitHsp {
                    q_seed: q,
                    s_seed: s,
                    q_start: u.q_start,
                    s_start: u.s_start,
                    length: u.length,
                    score: u.score,
                });
            }
            if options.window_size == 0 {
                entry.set_last_hit((end - (word - 1)).wrapping_add(diagonals.offset));
            } else if extended {
                entry.set_flag(1);
                entry.set_last_hit((end - (word - 1)).wrapping_add(diagonals.offset));
            } else {
                entry.set_last_hit(s.wrapping_add(diagonals.offset));
            }
        }
    };
    scan(&mut each);
    diagonals.finish_subject(subject.len().saturating_sub(1))?;
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
// ```c++
// static Int4 s_BlastDiagClear(BLAST_DiagTable * diag)
// {
//     Int4 i, n;
//     DiagStruct *diag_struct_array;
//
//     if (diag == NULL)
//         return 0;
//
//     n = diag->diag_array_length;
//
//     diag->offset = diag->window;
//
//     diag_struct_array = diag->hit_level_array;
//
//     for (i = 0; i < n; i++) {
//         diag_struct_array[i].flag = 0;
//         diag_struct_array[i].last_hit = -diag->window;
//         if (diag->hit_len_array) diag->hit_len_array[i] = 0;
//     }
//     return 0;
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
// ```c++
//     if (ewp->diag_table) {
//         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//             ewp->diag_table->offset = ewp->diag_table->window;
//             s_BlastDiagClear(ewp->diag_table);
//         } else {
//             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//         }
//     } else if (ewp->hash_table) {
// ```
#[cfg(test)]
mod tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    #[test]
    fn diagonal_overflow_threshold_clears_flags_and_negative_window() {
        let mut d = Diagonals::new(500, 40).unwrap();
        d.offset = i32::MAX / 4 - 1;
        let mut entry = DiagStruct::default();
        entry.set_last_hit(123);
        entry.set_flag(1);
        d.entries[0] = entry;
        d.finish_subject(83).unwrap();
        assert_eq!(d.offset, i32::MAX / 4 - 1 + 83 + 40);
        assert_eq!(d.entries[0], entry);
        d.finish_subject(83).unwrap();
        assert_eq!(d.offset, 40);
        assert!(d
            .entries
            .iter()
            .all(|x| x.last_hit() == -40 && x.flag() == 0));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:49-59
    // ```c++
    //
    //         if (diag_table)
    //         {
    //                 diag_array_length = 1;
    //                 /* What power of 2 is just longer than the query? */
    //                 while (diag_array_length < (qlen+window_size))
    //                 {
    //                         diag_array_length = diag_array_length << 1;
    //                 }
    //                 /* These are used in the word finders to shift and mask
    //                 rather than dividing and taking the remainder. */
    // ```
    #[test]
    fn diagonal_size_power_and_int4_overflow_boundary() {
        assert_eq!(Diagonals::new(512, 0).unwrap().entries.len(), 512);
        assert_eq!(Diagonals::new(512, 1).unwrap().entries.len(), 1024);
        assert!(Diagonals::new(i32::MAX as usize, 1).is_err());
        assert!(Diagonals::new(usize::MAX, 0).is_err());
    }
}

#[cfg(test)]
mod session_f_diagonal_boundary_tests {
    use super::{DiagStruct, Diagonals};

    /// An entry with `last_hit` and the flag set.
    pub(super) fn flagged(last_hit: i32) -> DiagStruct {
        let mut entry = DiagStruct::default();
        entry.set_last_hit(last_hit);
        entry.set_flag(1);
        entry
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
    // ```c++
    //         Int4 diag_array_length;
    //
    //         diag_table= (BLAST_DiagTable*) calloc(1, sizeof(BLAST_DiagTable));
    //
    //         if (diag_table)
    //         {
    //                 diag_array_length = 1;
    //                 /* What power of 2 is just longer than the query? */
    //                 while (diag_array_length < (qlen+window_size))
    //                 {
    //                         diag_array_length = diag_array_length << 1;
    //                 }
    //                 /* These are used in the word finders to shift and mask
    //                 rather than dividing and taking the remainder. */
    //                 diag_table->diag_array_length = diag_array_length;
    //                 diag_table->diag_mask = diag_array_length-1;
    //                 diag_table->multiple_hits = multiple_hits;
    //                 diag_table->offset = window_size;
    //                 diag_table->window = window_size;
    //         }
    // ```
    #[test]
    fn allocation_power_of_two_and_window_boundaries() {
        for (query, window, length) in [
            (63, 0, 64),
            (64, 0, 64),
            (65, 0, 128),
            (24, 40, 64),
            (25, 40, 128),
        ] {
            let d = Diagonals::new(query, window).unwrap();
            assert_eq!(d.entries.len(), length);
            assert_eq!(d.offset, window);
        }
        assert!(Diagonals::new(1, -1).is_err());
    }

    // No NCBI counterpart (buffer reuse, LOSAT_X_BXDIAGKEEP): a used table of another length or
    // window, made again by `x_remake`, is the table `new` makes.
    #[cfg(feature = "parallel")]
    #[test]
    fn a_remade_table_is_the_new_one() {
        for (from, to) in [
            ((100, 40), (30, 40)),
            ((30, 40), (100, 40)),
            ((70, 40), (70, 8)),
        ] {
            let mut table = Diagonals::new(from.0, from.1).unwrap();
            table.entries.fill(flagged(12_345));
            table.offset = 99;
            assert!(!table.x_fits(to.0, to.1));
            table.x_remake(to.0, to.1).unwrap();
            let new = Diagonals::new(to.0, to.1).unwrap();
            assert!(table.x_fits(to.0, to.1));
            assert_eq!(table.entries, new.entries);
            assert_eq!((table.offset, table.window), (new.offset, new.window));
        }
        assert!(Diagonals::new(8, 8).unwrap().x_remake(1, -1).is_err());
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    #[test]
    fn reset_occurs_at_entry_threshold_and_clears_stale_flags() {
        for window in [0, 1, 39, 40, 41] {
            let mut d = Diagonals::new(65, window).unwrap();
            d.entries.fill(flagged(123));
            d.offset = 536_870_910; // INT4_MAX / 4 - 1, independent literal boundary.
            d.finish_subject(1).unwrap();
            assert_eq!(d.offset, 536_870_911 + window);
            assert!(d.entries.iter().all(|x| *x == flagged(123)));
            d.finish_subject(5_000_000).unwrap();
            assert_eq!(d.offset, window);
            assert!(d.entries.iter().all(|x| *x == DiagStruct::clear(window)));
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    // Independent worker slots advance only for their assigned subjects, including
    // a partial final wave; no query state is shared between slots or batches.
    #[test]
    fn uneven_wave_slots_and_new_batch_keep_separate_offsets() {
        let mut slots = (0..4)
            .map(|_| Diagonals::new(128, 40).unwrap())
            .collect::<Vec<_>>();
        for (oid, length) in [10, 20, 30, 40, 50, 60, 70, 80, 90].into_iter().enumerate() {
            slots[oid % 4].finish_subject(length).unwrap();
        }
        assert_eq!(
            slots.iter().map(|x| x.offset).collect::<Vec<_>>(),
            [310, 200, 220, 240]
        );
        assert_eq!(Diagonals::new(128, 40).unwrap().offset, 40);
    }
}

#[cfg(test)]
mod session_f_word_finder_reset_tests {
    use super::session_f_diagonal_boundary_tests::flagged;
    use super::*;
    use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    // Expected SEED/INIT rows are frozen fresh pinned-NCBI trace bytes, not Rust output.
    #[test]
    fn fresh_retained_and_threshold_reset_word_finder_match_ncbi() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
        let q = parse_fasta(&std::fs::read(dir.join("query.fna")).unwrap(), false, false).unwrap();
        let s = parse_fasta(
            &std::fs::read(dir.join("subject.faa")).unwrap(),
            true,
            false,
        )
        .unwrap();
        let mut subject = s[0]
            .sequence
            .iter()
            .copied()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect::<Vec<_>>();
        subject.push(0);
        for window in [0, 1, 39, 40, 41] {
            let args = [
                "losat",
                "blastx",
                "-query",
                "query",
                "-subject",
                "subject",
                "-seg",
                "no",
                "-comp_based_stats",
                "0",
                "-window_size",
                &window.to_string(),
            ];
            let crate::cli::Commands::Blastx(args) =
                crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
                    .unwrap()
                    .command
            else {
                unreachable!()
            };
            let options = args.resolve().unwrap();
            let batch = prepare_queries(&q, &options).unwrap();
            let mut params = parameters::score_block(&batch);
            parameters::effective_lengths(
                &batch,
                &mut params,
                &options,
                s[0].sequence.len() as i64,
                1,
                None,
            )
            .unwrap();
            let lookup = Lookup::new(&batch, &params, &options).unwrap();
            parameters::subject_parameters(
                &batch,
                &mut params,
                &options,
                s[0].sequence.len() as i64,
                1,
                s[0].sequence.len(),
                None,
            )
            .unwrap();
            let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
            let expected =
                std::fs::read_to_string(dir.join(format!("window{window}.stage"))).unwrap();
            let seeds = expected
                .lines()
                .filter(|x| x.starts_with("SEED\t"))
                .map(|x| {
                    let f = x.split('\t').collect::<Vec<_>>();
                    (f[2].parse::<i32>().unwrap(), f[3].parse::<i32>().unwrap())
                })
                .collect::<Vec<_>>();
            let initial = expected
                .lines()
                .filter(|x| x.starts_with("INIT\t"))
                .map(|x| {
                    let f = x.split('\t').collect::<Vec<_>>();
                    InitHsp {
                        q_seed: f[3].parse().unwrap(),
                        s_seed: f[4].parse().unwrap(),
                        q_start: f[5].parse().unwrap(),
                        s_start: f[6].parse().unwrap(),
                        length: f[7].parse().unwrap(),
                        score: f[8].parse().unwrap(),
                    }
                })
                .collect::<Vec<_>>();
            let mut d = Diagonals::new(length, window).unwrap();
            for state in 0..4 {
                if state == 2 {
                    d.offset = 536_870_910;
                    d.entries.fill(flagged(123));
                }
                if state == 3 {
                    d.offset = 536_870_911;
                    d.entries.fill(flagged(536_870_900));
                    d.finish_subject(83).unwrap();
                }
                let mut actual = Vec::new();
                let hits = word_finder(&batch, &params, &options, &subject, &lookup, &mut d, |x| {
                    actual.extend_from_slice(x)
                })
                .unwrap();
                assert_eq!(actual, seeds, "window={window} state={state}");
                assert_eq!(hits, initial, "window={window} state={state}");
            }
        }
    }
}
