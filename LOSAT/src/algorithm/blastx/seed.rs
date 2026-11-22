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
    s_left_off: i32,
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
// The seed loop on a per-thread offset-pair buffer, read in place (a child module, so that it uses
// the diagonal table's fields).
#[path = "x_seed_scan.rs"]
pub mod x_seed_scan;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-511
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
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
// The entries are NCBI's bit-field (`DiagStruct`: `last_hit` truncated to 31 bits on write and
// sign-extended on read), so a `s + diag_offset` past 2^30 wraps as in NCBI. That happens at a
// `-window_size` near 2^29 from the second subject on (the offset grows by length + window per
// subject until INT4_MAX / 4); NCBI 2.17.0 then drops that subject's two-hit seeds, and so does
// LOSAT. The sums and differences on the entries wrap (two's complement i32), as NCBI's do.
pub struct Diagonals {
    entries: Vec<DiagStruct>,
    offset: i32,
    window: i32,
    /// `-db_soft_mask` and `-db_hard_mask`: the worker's reader of the stored masks
    /// (`search_db_masked.rs`), made at its first masked sequence; a search without them keeps
    /// `None` (one word in the table).
    pub(crate) x_db_masks: Option<Box<crate::algorithm::common::db_mask_search::SubjectMasks>>,
}
impl Diagonals {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
    pub fn new(query_length: usize, window: i32) -> Result<Self> {
        Ok(Self {
            entries: vec![DiagStruct::default(); Self::array_length(query_length, window)?],
            offset: window,
            window,
            x_db_masks: None,
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
        self.x_db_masks = None;
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
    /// [`Diagonals::finish_chunk`] after a subject of `length` letters.
    pub fn finish_subject(&mut self, length: usize) -> Result<()> {
        self.finish_chunk(i32::try_from(length)?)
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
    /// `Blast_ExtendWordExit` after a subject, or a hard chunk of a `-db_hard_mask` sequence
    /// (`search_db_masked.rs`), whose length is NCBI's `Int4` (`core::subject_split::HardChunk::
    /// length`: negative between overlapping masks): the table cleared where the offset reached
    /// INT4_MAX / 4, else the offset moved by the length and the window. The one owner of that
    /// update.
    pub(crate) fn finish_chunk(&mut self, length: i32) -> Result<()> {
        if self.offset >= i32::MAX / 4 {
            self.clear();
        } else {
            self.offset = self
                .offset
                .checked_add(length)
                .and_then(|v| v.checked_add(self.window))
                .ok_or_else(|| anyhow::anyhow!("BLASTX diagonal offset exceeds Int4"))?;
        }
        Ok(())
    }
    // No NCBI counterpart: PD-LOSAT-NCBI-DEFECTS exception 23 (the Owner's answer Q5 = B of M4b).
    /// The offset moved up by `by` between two sequences of a `-db_hard_mask` search
    /// (`search_db_masked.rs`), so that no entry of one reaches the other's hits
    /// (`common::hard_chunk_walk::DiagonalWalk::start`, `DiagonalWalk::after`: below INT4_MAX / 4
    /// after the move, else they clear the table instead).
    pub(crate) fn move_offset(&mut self, by: i32) -> Result<()> {
        self.offset = self
            .offset
            .checked_add(by)
            .ok_or_else(|| anyhow::anyhow!("BLASTX diagonal offset exceeds Int4"))?;
        Ok(())
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
    /// Clears the table and sets the offset to the window (`s_BlastDiagClear`): where the offset
    /// reached INT4_MAX / 4 (`finish_chunk`), and in a `-db_hard_mask` search between two
    /// sequences where moving the offset to keep them apart would take it there
    /// (`search_db_masked.rs`, `common::hard_chunk_walk::DiagonalWalk::start`,
    /// `DiagonalWalk::after`).
    pub(crate) fn clear(&mut self) {
        self.offset = self.window;
        self.entries.fill(DiagStruct::clear(self.window));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66, 135-150
    /// The table as NCBI's thread makes it (`BlastExtendWordNew`: every entry 0 without the flag,
    /// the offset the window), for the first searched sequence of a `-db_hard_mask` search
    /// (`search_db_masked.rs`), whose chunks can go back below its start onto those entries
    /// (`common::hard_chunk_walk::DiagonalWalk::start_on_new_table`).
    pub(crate) fn new_table(&mut self) {
        self.offset = self.window;
        self.entries.fill(DiagStruct::default());
    }
    /// The diagonal offset of the next subject.
    pub(crate) fn offset(&self) -> i32 {
        self.offset
    }
    /// The entries, the mask of a diagonal's index and the diagonal offset, which the two-hit test
    /// of a hard chunk's one word reads (`blastp::hard_chunk_scan::OneWord::test`).
    pub(crate) fn table(&mut self) -> (&mut [DiagStruct], i32, i32) {
        let (mask, offset) = (self.entries.len() as i32 - 1, self.offset);
        (self.entries.as_mut_slice(), mask, offset)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:273-313
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
    /// [`Lookup::x_new_for`] with the cutoff of the search's options (`x_seed_diag::x_narrow_cutoff`).
    pub(crate) fn x_new_in_parts(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
        parts: usize,
    ) -> Result<Self> {
        Lookup::x_new_for(
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
/// EXPERIMENT (LOSAT_X_BPCLUT): the compressed table of a query batch (`-word_size` 5 to 7), flat as
/// TBLASTN's with `u16` query offsets (`i32` past `X_U16_LETTERS` letters;
/// `x_bx_build_compressed_lookup`, on vectors kept between builds with `LOSAT_X_BXCLUTPOOL`), on one
/// thread, or in `parts` above 1 (`LOSAT_X_BXBUILDPAR`, the same table). One
/// table for the batch, so its one
/// range is all the segments. Words that fill NCBI's overflow banks are refused, as BLASTP's
/// (NCBI writes past its banks). Cold and not inlined, so that the standard table's build keeps
/// its code.
#[cold]
#[inline(never)]
fn x_compressed_lookup(
    query: &[u8],
    segments: &[(i32, i32)],
    threshold: f64,
    word: i32,
    parts: usize,
) -> Result<Lookup> {
    // The build of the word length, as a constant of its instance (`x_bx_build_compressed_lookup`).
    let scaled = x_compressed_scaled_threshold(threshold);
    let tables = [0..segments.len()];
    match word {
        5 => x_bx_build_compressed_lookup::<5>(scaled, query, segments, &tables, parts),
        6 => x_bx_build_compressed_lookup::<6>(scaled, query, segments, &tables, parts),
        7 => x_bx_build_compressed_lookup::<7>(scaled, query, segments, &tables, parts),
        _ => None,
    }
    .ok_or_else(|| anyhow::anyhow!("unsupported BLASTX compressed lookup"))?
    .map(Lookup::Compressed)
    .map_err(|()| {
        anyhow::anyhow!(
            "-threshold {threshold} with the compressed lookup table of -word_size {word} is not supported by LOSAT's BLASTX: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap"
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
        .map(|(c, p)| {
            let mut aa = vec![0];
            aa.extend_from_slice(&query[c.offset..c.offset + c.length]);
            aa.push(0);
            QueryContext {
                q_idx: c.query_index as u32,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1349-1351
    // Only the tests build a table this way; the searches use `x_new_for` (`BxLookup::x_new_in_parts`).
    #[cfg(test)]
    pub(crate) fn new(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
    ) -> Result<Self> {
        match Self::x_new_for(batch, parameters, options, 1, None)? {
            BxLookup::Wide(lookup) => Ok(lookup),
            BxLookup::Narrow(_) => unreachable!("a build without a cutoff takes Int4 hits"),
        }
    }

    /// The table of a query batch: the search's matrix and word size choose the builds once per
    /// query batch. BLOSUM62's of `-word_size 3` and the compressed table's (5 to 7) are
    /// `x_new_as`'s; another matrix's and `-word_size` 2 and 4 of the standard table are
    /// `x_new_of`'s (the same tables with the matrix's scores and the word length).
    fn x_new_for(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
        parts: usize,
        narrow_cutoff: Option<u64>,
    ) -> Result<BxLookup> {
        if options.scoring != ScoringMatrix::Blosum62
            || (!options.compressed_lookup() && options.word_size != 3)
        {
            return Self::x_new_of(batch, parameters, options, parts, narrow_cutoff);
        }
        Self::x_new_as(batch, parameters, options, parts, narrow_cutoff)
    }

    /// EXPERIMENT (LOSAT_X_BXBUILDPAR): the table of `Lookup::new` (the tests' one-task build)
    /// with a pooled standard table's passes on `parts` tasks of the pool it runs on
    /// (`tblastx/lookup/x_lut_direct.rs` `x_partitioned_passes`); the table is the same. The
    /// standard table takes 16-bit hits (a [`NarrowLookup`]) where the lookup query fits and its
    /// hits per used cell, in hundredths, reach `narrow_cutoff` (`None`: never).
    pub(crate) fn x_new_as(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
        parts: usize,
        narrow_cutoff: Option<u64>,
    ) -> Result<BxLookup> {
        let query = batch.sequence_start[1..].to_vec();
        if options.compressed_lookup() {
            return x_compressed_lookup(
                &query,
                &batch.lookup_segments,
                options.threshold,
                options.word_size,
                parts,
            )
            .map(BxLookup::Wide);
        }
        ensure!(
            options.word_size == 3,
            "unsupported BLASTX lookup word size"
        );
        Self::x_standard_with(
            batch,
            parameters,
            options,
            query,
            true,
            |query, segments, contexts, threshold, pooled| {
                build_ncbi_lookup_from_prepared(
                    query,
                    segments,
                    contexts,
                    threshold,
                    pooled,
                    parts,
                    narrow_cutoff,
                )
            },
        )
    }

    /// The standard table of [`Lookup::x_new_as`] and [`Lookup::x_new_of`], with `build` the
    /// matrix's builder (its query, lookup segments, contexts, threshold and whether it takes the
    /// pool). `poolable`: whether the pooled build covers the table's backbone index (16 bits:
    /// words of 2 and 3 letters). A table of 4 letters (the joint builder, item SC9 of line NP)
    /// takes no kept vectors, so it gives none back: no build of the search would take them.
    #[inline(always)]
    fn x_standard_with(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
        query: Vec<u8>,
        poolable: bool,
        build: impl FnOnce(Vec<u8>, Vec<(i32, i32)>, Vec<QueryContext>, i32, bool) -> XBxTable,
    ) -> Result<BxLookup> {
        let contexts = lookup_contexts(batch, parameters, &query);
        // Dispatch point (LOSAT_X_BXPOOLN1, in the strict set): a search on one thread builds its
        // tables without the pool (tblastx/lookup/x_lut_direct.rs). No NCBI counterpart.
        let pooled = poolable
            && crate::algorithm::tblastx::lookup::x_bx_pooled_for(
                usize::try_from(options.num_threads).unwrap_or(1),
            );
        let kept = x_seed_scan::XKeptPresence::x_built(pooled);
        Ok(
            match build(
                query,
                batch.lookup_segments.clone(),
                contexts,
                // `lookup->threshold` is `opt->threshold` cast to `Int4` (blast_aalookup.c:241), converted as
                // the x86-64 build does (`INT_MIN` past `Int4`).
                crate::core::blast_util::ncbi_int4_from_double(options.threshold),
                pooled,
            ) {
                XBxTable::Wide(table) => BxLookup::Wide(Self::Standard(table, kept)),
                XBxTable::Narrow(table) => BxLookup::Narrow(NarrowLookup { table, kept }),
            },
        )
    }
    /// [`Lookup::x_new_as`] with another matrix than BLOSUM62 (items SC6 and SC6b of line NP) or the
    /// standard table of another word size than 3 (`-word_size` 2 and 4, item SC9): the compressed
    /// table of `-word_size` 5 to 7 with the matrix's scores (`x_bx_build_compressed_lookup_of`:
    /// BLOSUM62's builders with the matrix's words, `Int4` hits, in `parts`), the
    /// standard one as BLOSUM62's
    /// (`build_ncbi_lookup_from_prepared_of`: the pool, `parts`, 16-bit hits from `narrow_cutoff`)
    /// with its scores and word length. Not inlined, so that BLOSUM62's builds keep their code.
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
            x_bx_build_compressed_lookup_of, XCompressedLookup,
        };
        let matrix = options.scoring;
        let query = batch.sequence_start[1..].to_vec();
        let word = usize::try_from(options.word_size)?;
        if options.compressed_lookup() {
            let segments = &batch.lookup_segments;
            return x_bx_build_compressed_lookup_of(
                matrix,
                word,
                x_compressed_scaled_threshold(options.threshold),
                &query,
                segments,
                &[0..segments.len()],
                parts,
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
                    "-threshold {} with the compressed lookup table of -word_size {word} is not supported by LOSAT's BLASTX: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap",
                    options.threshold
                )
            });
        }
        ensure!(
            (2..=4).contains(&word),
            "unsupported BLASTX lookup word size"
        );
        Self::x_standard_with(
            batch,
            parameters,
            options,
            query,
            crate::algorithm::tblastx::lookup::x_pair_lut::x_direct_covers(word),
            |query, segments, contexts, threshold, pooled| {
                crate::algorithm::tblastx::lookup::build_ncbi_lookup_from_prepared_of(
                    matrix,
                    word,
                    query,
                    segments,
                    contexts,
                    threshold,
                    pooled,
                    parts,
                    narrow_cutoff,
                )
            },
        )
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
    pub fn capacity(&self) -> usize {
        4096 + match self {
            Self::Standard(l, _) => l.longest_chain as usize,
            Self::Compressed(XBxCompressedLookup::Reference(l)) => l.longest_chain as usize,
            Self::Compressed(XBxCompressedLookup::Flat(l)) => l.longest_chain as usize,
            Self::Compressed(XBxCompressedLookup::Wide(l)) => l.longest_chain as usize,
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-769
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:499-500
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:774-800
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-515
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
// word_finder is the original port of this loop for one subject; LOSAT_X_BXDIAGSHADOW and the
// unit tests of the faster loops compare them with it.
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
// Only the tests call it: the reference that the tests of `x_seed_diag` compare the searches'
// word finders with.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(super) fn word_finder_bx_of<S: BxScoring>(
    scoring: S,
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
// Only the tests call it: the reference of the tests of `x_seed_diag` and `x_soft_subject`.
#[cfg(test)]
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
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-510, 566-573
                // The scan's query offset is unsigned: the context test subtracts the distance in
                // unsigned arithmetic, and a distance past the query's start (an entry of an
                // earlier hard chunk of the sequence, series NP, M4b-3b) passes it.
                if (q.wrapping_sub(diff) as u32) < batch.contexts[context].offset as u32 {
                    entry.set_last_hit(s.wrapping_add(diagonals.offset));
                    continue;
                }
                // The first hit's end is signed: negative where it lies before a hard chunk's
                // start.
                let Some(result) = bx_two_hit(
                    scoring,
                    query,
                    subject,
                    previous + word,
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
    sort_initial_hits(&mut hits);
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234-235
/// The order of a subject's initial hits for its gapped stage (by score, then the subject's start,
/// the length and the query's start): the reference word finder's (`word_finder`) and a hard
/// chunk's one word's (`search_db_masked.rs`).
pub(crate) fn sort_initial_hits(hits: &mut [InitHsp]) {
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
#[cfg(test)]
#[path = "seed_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "seed_session_f_diagonal_boundary_tests.rs"]
mod session_f_diagonal_boundary_tests;

#[cfg(test)]
#[path = "seed_session_f_word_finder_reset_tests.rs"]
mod session_f_word_finder_reset_tests;
