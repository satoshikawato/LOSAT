//! Internal NCBI TBLASTN two-hit WordFinder of the (frame, chunk) units of the preliminary search
//! (`XTwoHit`).

#[cfg(test)]
use super::search_seed::encode_tblastn_lookup_query;
use crate::algorithm::blastp::extension::{
    extend_two_hit, extend_two_hit_outside, BlastpOneHitUngappedResult, BlastpTwoHitUngappedResult,
    ChunkLetters,
};
use crate::algorithm::blastp::x_ungapped_fast;
use crate::config::ScoringMatrix;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:142-163
// Mirrors `BlastUngappedData` (query start, subject start, length, and score) and
// `BlastInitHSP` (seed offset pair plus ungapped extension data).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct InitHsp {
    pub frame: i8,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:241-250,572-586:
    // backup.offset is retained until HSP offsets are adjusted before merge.
    pub chunk_offset: u32,
    pub q_seed: u32,
    pub s_seed: u32,
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65,162-175
pub(super) struct Diagonals {
    entries: XDiagCells,
    mask: i32,
    offset: i32,
    window: i32,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
/// The cells of a diagonal table: NCBI's `DiagStruct` in four bytes (`LOSAT_X_TBNDIAGPACK`: last_hit
/// in bits 0-30, flag in bit 31, `x_tbn_hitbf::XDiagCell`), or the eight of `(last_hit, flag)`.
#[derive(Debug, PartialEq)]
pub(super) enum XDiagCells {
    Packed(Vec<u32>),
    Pair(Vec<(i32, bool)>),
}

impl XDiagCells {
    fn len(&self) -> usize {
        match self {
            Self::Packed(cells) => cells.len(),
            Self::Pair(cells) => cells.len(),
        }
    }

    fn packed(&self) -> bool {
        matches!(self, Self::Packed(_))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
    // `s_BlastDiagClear` sets every cell's `last_hit` to `-window` and clears `flag`. A packed
    // table never runs a lead subject (its offsets and hits are non-negative), and for hits at
    // subject offset 0 or more a cleared cell of 0 and one of `-window` give the same outcome:
    // both lie at least a full window behind every hit when `offset == window`.
    fn clear(&mut self, window: i32) {
        match self {
            Self::Packed(cells) => cells.fill(0),
            Self::Pair(cells) => cells.fill((-window, false)),
        }
    }
}

// No NCBI counterpart: reads the LOSAT_X_TBNDIAGPACK variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNDIAGPACK): a switch of the strict set (`utils/xswitch.rs`), on unless it
/// is `0`: the diagonal tables of the two-hit word finder hold NCBI's four-byte `DiagStruct`
/// (`XDiagCells::Packed`) where the table's values fit it (`x_packs`).
fn x_tbn_diag_pack() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNDIAGPACK"))
}

/// No NCBI counterpart: the largest window of a packed table (`x_packs`), 2^28.
const X_PACK_MAX_WINDOW: i32 = 1 << 28;

// No NCBI counterpart: chooses the layout of a diagonal table; both hold the same values.
/// Whether a table of `window` and `word_size` is packed (`LOSAT_X_TBNDIAGPACK`). Its values then
/// fit 31 bits without a sign: `last_hit` starts at 0 and is set to `s + offset` or to
/// `s_last_off - (wordsize - 1) + offset` (`XTwoHit::pairs`, `s <= s_last_off <=` the unit's
/// length), with `offset` from `window` (at least `wordsize + 1`, `x_tbn_hitbf::fits`) to below
/// INT4_MAX / 4 + window + the longest unit (`Diagonals::finish_frame`), the units at most
/// `MAX_DBSEQ_LEN` (5,000,000) residues long (`search_seed.rs`) and the window at most
/// `X_PACK_MAX_WINDOW`: below 8.3 x 10^8 < 2^30, also inside the range of NCBI's signed
/// `last_hit:31` (-2^30 to 2^30 - 1), so both tables hold NCBI's values. Bit 31 is free for the
/// flag, and a flag-0 cell read as `i32` is its `last_hit` (`x_tbn_hitbf::x_hits_bf`). TBLASTN
/// takes window 40 with word size 3 and 60 with 2 (`args.rs` `search_settings`), both packed.
fn x_packs(window: i32, word_size: usize) -> bool {
    x_tbn_diag_pack() && super::x_tbn_hitbf::fits(window, word_size) && window <= X_PACK_MAX_WINDOW
}

impl Diagonals {
    fn new(query_length: usize, window: i32, packed: bool) -> Self {
        let size = (query_length + window as usize).next_power_of_two();
        Self {
            entries: if packed {
                XDiagCells::Packed(vec![0; size])
            } else {
                XDiagCells::Pair(vec![(0, false); size])
            },
            mask: size as i32 - 1,
            offset: window,
            window,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
    // `subject_length` is NCBI's Int4 chunk length (negative between overlapping hard masks).
    fn finish_frame(&mut self, subject_length: i32) {
        if self.offset >= i32::MAX / 4 {
            self.offset = self.window;
            self.entries.clear(self.window);
        } else {
            self.offset = self
                .offset
                .wrapping_add(subject_length)
                .wrapping_add(self.window);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66,87-107,131-144
// Sets a lead subject's diagonal table (`Pair` cells) to NCBI's one-thread start offset and
// unwritten-cell value (`TableStart::New`, series NP, M4b-T3).
impl crate::algorithm::common::hard_chunk_walk::HardTable for Diagonals {
    type Error = anyhow::Error;

    fn ncbi_table(
        &mut self,
        start: crate::algorithm::common::hard_chunk_walk::NcbiStart,
        window: i32,
    ) {
        debug_assert_eq!(window, self.window);
        self.offset = start.offset;
        match &mut self.entries {
            XDiagCells::Pair(cells) => cells.fill((start.unwritten(self.window), false)),
            XDiagCells::Packed(_) => {
                debug_assert!(
                    false,
                    "packed diagonal table never runs a soft-mask lead subject"
                );
            }
        }
    }

    fn clear_table(&mut self, window: i32) {
        debug_assert_eq!(window, self.window);
        self.offset = self.window;
        self.entries.clear(self.window);
    }

    fn move_offset(&mut self, by: i32) -> anyhow::Result<()> {
        self.offset = self
            .offset
            .checked_add(by)
            .ok_or_else(|| anyhow::anyhow!("TBLASTN diagonal offset exceeds Int4"))?;
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-518, 776-779
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
/// TBLASTN's `(last_hit, flag)` diagonal ring as a [`DiagTable`] for
/// `common::hard_chunk_scan::OneWord`.
pub(super) struct TbnPairRing<'a> {
    entries: &'a mut [(i32, bool)],
    mask: i32,
}

impl crate::algorithm::common::diag_entry::DiagTable for TbnPairRing<'_> {
    type Cell = (i32, bool);

    #[inline(always)]
    fn two_hit_entry(&mut self, q: i32, s: i32) -> (&mut (i32, bool), i32) {
        (
            &mut self.entries[(q.wrapping_sub(s) & self.mask) as usize],
            0,
        )
    }

    #[inline(always)]
    fn one_hit_entry(&mut self, q: i32, s: i32) -> &mut (i32, bool) {
        &mut self.entries[(s.wrapping_sub(q) & self.mask) as usize]
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-547
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
// XTwoHit is the two-hit pass of NCBI's word finder with its own diagonal array per (frame, chunk)
// unit. Each unit starts from the state of a fresh table, which is what the offset update at the
// end of the previous WordFinder call gives (see the comment above), so units can run in any order
// or at the same time.
// The two-hit pass of one (frame, chunk) unit on its
// own diagonal array.  NCBI c++/src/algo/blast/core/blast_extend.c:162-175 ends a
// WordFinder call by adding `subject_length` and `window` to `diag_table->offset`, which
// makes every entry of the previous call older than a full window: the first hit
// of the next call on that diagonal only records `last_hit` (aa_ungapped.c:
// 516-547), exactly as on a cleared array.  A unit therefore depends on nothing
// that an earlier unit left behind.
pub(super) struct XTwoHit<'a> {
    diagonals: Diagonals,
    pub query_sequence: &'a [u8],
    pub context_offsets: &'a [i32],
    pub x_dropoffs: &'a [i32],
    pub cutoff_scores: &'a [i32],
    pub matrix: ScoringMatrix,
    pub window: i32,
    pub word_size: usize,
    /// EXPERIMENT (LOSAT_X_TBNUNGAP): the extension takes BLASTP's byte-table walks
    /// ([`x_tbn_ungap`]).
    x_ungap_fast: bool,
    /// The table of those walks for a matrix other than BLOSUM62 (`x_tbn_ungap`); `None` for
    /// BLOSUM62 and for the reference.
    x_ungap_matrix: Option<&'static x_ungapped_fast::MatrixTable>,
    /// The letters around the unit's chunk for the extensions of its hits outside the chunk
    /// (`x_set_around`); `None` for a unit whose scan gives no such hit.
    x_around: Option<XAround>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:300-306 (cited at
// `ChunkLetters`)
/// A copy of the letters around a unit's chunk (`ChunkLetters`): made only for a unit whose scan
/// reads letters before its chunk (masks stored out of order), so that the finder, which lives
/// longer than a unit's letters, holds them for that unit.
struct XAround {
    letters: Vec<u8>,
    lead: usize,
    length: usize,
}

impl XAround {
    fn letters(&self) -> ChunkLetters<'_> {
        ChunkLetters {
            letters: &self.letters,
            lead: self.lead,
            length: self.length,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1088-1096
// EXPERIMENT (LOSAT_X_TBNUNGAP): the two-hit extension of tblastn's word finder (`XTwoHit::pairs`)
// with BLASTP's byte-table walks (`blastp/x_ungapped_fast.rs`, whose module doc
// gives why the result is the reference's): the same statements in the same order, a residue pair
// scored from one signed byte (`q << 5 | s`) instead of the matrix row, and the bounds of each walk
// computed from the slice lengths as the reference computes them, so the walks stop where the
// reference's do on every unit, also on a chunk that ends inside a frame. An extension whose codes
// or X-drop the byte table does not cover (a code >= 32, an X-drop outside 1..=127, a word longer
// than the table takes: 11 residues for BLOSUM62, `MatrixTable::max_word_size` for another
// matrix) runs the reference. BLOSUM62, the matrix of every default tblastn search,
// takes `two_hit` (its table, where `build_table` gives one: every BLOSUM62 score fits a byte;
// were it not so, the walks hand every call to the reference). Another matrix takes `two_hit_of`
// with its own table (`matrix_table`, whose module comment gives the longest word the table
// takes; a matrix whose scores do not fit a byte keeps the reference). A strict-set switch
// (`utils/xswitch.rs`): on unless `LOSAT_X_TBNUNGAP=0`.
/// Whether a word finder with `matrix` takes BLASTP's byte-table walks (switch on), and the table
/// of a matrix other than BLOSUM62 (`None` for BLOSUM62's own walks).
fn x_tbn_ungap(matrix: ScoringMatrix) -> (bool, Option<&'static x_ungapped_fast::MatrixTable>) {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if !*ON.get_or_init(|| {
        crate::utils::xswitch::on("LOSAT_X_TBNUNGAP") && x_ungapped_fast::build_table().is_some()
    }) {
        return (false, None);
    }
    let table = (matrix != ScoringMatrix::Blosum62)
        .then(|| x_ungapped_fast::matrix_table(matrix))
        .flatten();
    let fast = matrix == ScoringMatrix::Blosum62 || table.is_some();
    crate::utils::xreach::x_reach!("LOSAT_X_TBNUNGAP", fast);
    (fast, table)
}

/// The byte-table walks of the two-hit loop (`XTwoHit::pairs_on`): BLOSUM62's (`XBlosum62Walks`)
/// or another matrix's (`&MatrixTable`), or the extension of a PSSM query's rows
/// (`score::XPssmWalks`, tblastn `-in_pssm`).
pub(super) trait XUngapWalks: Copy {
    /// `x_ungapped_fast::two_hit` with this table.
    #[allow(clippy::too_many_arguments)]
    fn two_hit(
        self,
        query: &[u8],
        subject: &[u8],
        s_left_off: i32,
        s_right_off: usize,
        q_right_off: usize,
        dropoff: i32,
        word_size: usize,
    ) -> Option<Option<crate::algorithm::blastp::extension::BlastpTwoHitUngappedResult>>;

    /// The extension of a hit outside its chunk (`XTwoHit::x_two_hit_around`) with these walks:
    /// a PSSM's own (`score::XPssmWalks`); `None` for a matrix's walks, whose extension there is
    /// the matrix's (`blastp::extension::extend_two_hit_outside`, the reference of the walks).
    #[allow(clippy::too_many_arguments)]
    fn two_hit_around(
        self,
        _chunk: ChunkLetters<'_>,
        _s_left_off: i32,
        _s_right_off: i32,
        _q_right_off: i32,
        _dropoff: i32,
        _word_size: i32,
    ) -> Option<BlastpTwoHitUngappedResult> {
        None
    }
}

/// BLOSUM62's walks (`x_ungapped_fast::two_hit`).
#[derive(Clone, Copy)]
struct XBlosum62Walks;

impl XUngapWalks for XBlosum62Walks {
    #[inline(always)]
    fn two_hit(
        self,
        query: &[u8],
        subject: &[u8],
        s_left_off: i32,
        s_right_off: usize,
        q_right_off: usize,
        dropoff: i32,
        word_size: usize,
    ) -> Option<Option<crate::algorithm::blastp::extension::BlastpTwoHitUngappedResult>> {
        x_ungapped_fast::two_hit(
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

impl XUngapWalks for &'static x_ungapped_fast::MatrixTable {
    #[inline(always)]
    fn two_hit(
        self,
        query: &[u8],
        subject: &[u8],
        s_left_off: i32,
        s_right_off: usize,
        q_right_off: usize,
        dropoff: i32,
        word_size: usize,
    ) -> Option<Option<crate::algorithm::blastp::extension::BlastpTwoHitUngappedResult>> {
        x_ungapped_fast::two_hit_of(
            self,
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

impl<'a> XTwoHit<'a> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1003,1377-1379,1412-1413
    // NCBI's search thread makes its diagonal table once, before its loop over the subjects, and
    // keeps it for every subject. A table of the same length, window and layout is taken over as
    // it is: `exit` after each unit leaves every entry older than a full window (see `exit`).
    /// `new` on `table` (a table an earlier unit left) when it has the length, window and layout
    /// (`x_packs`) this query needs, otherwise on a new table.
    pub(super) fn with_table(
        table: Option<Diagonals>,
        query_sequence: &'a [u8],
        context_offsets: &'a [i32],
        x_dropoffs: &'a [i32],
        cutoff_scores: &'a [i32],
        matrix: ScoringMatrix,
        window: i32,
        word_size: usize,
    ) -> Self {
        let packed = x_packs(window, word_size);
        Self::x_with_layout(
            packed,
            table,
            query_sequence,
            context_offsets,
            x_dropoffs,
            cutoff_scores,
            matrix,
            window,
            word_size,
        )
    }

    /// `with_table` with a table packed or not (`x_packs`; a test may take either).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn x_with_layout(
        packed: bool,
        table: Option<Diagonals>,
        query_sequence: &'a [u8],
        context_offsets: &'a [i32],
        x_dropoffs: &'a [i32],
        cutoff_scores: &'a [i32],
        matrix: ScoringMatrix,
        window: i32,
        word_size: usize,
    ) -> Self {
        let query_length = query_sequence.len().saturating_sub(1);
        let size = (query_length + window as usize).next_power_of_two();
        let diagonals = match table {
            Some(table)
                if table.entries.len() == size
                    && table.window == window
                    && table.entries.packed() == packed =>
            {
                table
            }
            _ => Diagonals::new(query_length, window, packed),
        };
        let (x_ungap_fast, x_ungap_matrix) = x_tbn_ungap(matrix);
        Self {
            diagonals,
            query_sequence,
            context_offsets,
            x_dropoffs,
            cutoff_scores,
            matrix,
            window,
            word_size,
            x_ungap_fast,
            x_ungap_matrix,
            x_around: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:300-306
    /// The letters around the next unit's chunk (`ChunkLetters`) for the extensions of its hits
    /// outside the chunk (`x_two_hit_around`, `x_one_hit_around`), or `None` for a unit whose
    /// scan gives no such hit (a copy; `None` frees it).
    pub(super) fn x_set_around(&mut self, around: Option<ChunkLetters<'_>>) {
        self.x_around = around.map(|around| XAround {
            letters: around.letters.to_vec(),
            lead: around.lead,
            length: around.length,
        });
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-583,1089-1158 (cited
    // at `blastp::extension::extend_two_hit_outside`)
    /// The extension of `pairs_on` for a hit whose first hit (`s_left_off`) lies before its
    /// unit's chunk, on the letters around the chunk (`x_set_around`): the walks' own (a PSSM's)
    /// or the matrix's. `None` for a unit without those letters (none has such a hit).
    #[cold]
    #[inline(never)]
    fn x_two_hit_around<W: XUngapWalks>(
        &self,
        walks: W,
        s_left_off: i32,
        s_right_off: i32,
        q_right_off: i32,
        dropoff: i32,
    ) -> Option<BlastpTwoHitUngappedResult> {
        let around = self.x_around.as_ref()?.letters();
        let word_size = self.word_size as i32;
        let walked = if self.x_ungap_fast {
            walks.two_hit_around(
                around,
                s_left_off,
                s_right_off,
                q_right_off,
                dropoff,
                word_size,
            )
        } else {
            None
        };
        Some(walked.unwrap_or_else(|| {
            extend_two_hit_outside(
                self.matrix,
                self.query_sequence,
                around,
                s_left_off,
                s_right_off,
                q_right_off,
                dropoff,
                word_size,
            )
        }))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:785-791,1020-1083 (cited
    // at `blastp::extension::extend_one_hit_outside`)
    /// The extension of `one_hit_on` for a word (at `s_off`) before its unit's chunk, on the
    /// letters around the chunk (`x_set_around`); `None` for a unit without them.
    #[cold]
    #[inline(never)]
    fn x_one_hit_around<S: super::score::TbnScore>(
        &self,
        score: S,
        q_off: i32,
        s_off: i32,
        dropoff: i32,
    ) -> Option<BlastpOneHitUngappedResult> {
        let around = self.x_around.as_ref()?.letters();
        score.one_hit_around(
            self.query_sequence,
            around,
            q_off,
            s_off,
            dropoff,
            self.word_size as i32,
            self.matrix,
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:614,807
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:166-172
    // The end of NCBI's WordFinder call for one (frame, chunk). Every entry the unit wrote holds
    // at most `offset + subject_length - 1`, so after the step each one lies more than a window
    // behind any hit of the next unit: its first hit on a diagonal only records `last_hit` and
    // clears `flag` (the tests of `pairs`), as on a table made afresh by `Diagonals::new`. The
    // next unit sees this instead of a new table.
    /// End a unit of `subject_length` residues (NCBI's Int4, negative between overlapping hard
    /// masks) as NCBI does (the table is kept).
    pub(super) fn exit(&mut self, subject_length: i32) {
        self.diagonals.finish_frame(subject_length);
    }

    /// Current diagonal offset of the table (NCBI `diag_table->offset`).
    pub(super) fn offset(&self) -> i32 {
        self.diagonals.offset
    }

    /// Borrows the `(last_hit, flag)` diagonal ring with its current offset for `OneWord` (`None`
    /// on a packed table).
    pub(super) fn pair_ring(&mut self) -> Option<(TbnPairRing<'_>, i32)> {
        let mask = self.diagonals.mask;
        let offset = self.diagonals.offset;
        match &mut self.diagonals.entries {
            XDiagCells::Pair(entries) => Some((TbnPairRing { entries, mask }, offset)),
            XDiagCells::Packed(_) => None,
        }
    }

    /// Applies `start` to this word finder's diagonal table (`HardTable`, series NP, M4b-T3).
    pub(super) fn apply_table_start(
        &mut self,
        start: crate::algorithm::common::hard_chunk_walk::TableStart,
    ) -> anyhow::Result<()> {
        start.apply(&mut self.diagonals, self.window)
    }

    /// No NCBI counterpart: hands the diagonal table back for `with_table`; it does not change
    /// any value NCBI computes.
    pub(super) fn into_table(self) -> Diagonals {
        self.diagonals
    }

    /// The diagonal table as `(last_hit, flag)`, for tests.
    #[cfg(test)]
    pub(super) fn x_table(&self) -> Vec<(i32, bool)> {
        use super::x_tbn_hitbf::XDiagCell;
        match &self.diagonals.entries {
            XDiagCells::Packed(cells) => cells.iter().map(|cell| cell.get()).collect(),
            XDiagCells::Pair(cells) => cells.clone(),
        }
    }

    /// The number of cells of the diagonal table, for tests.
    #[cfg(test)]
    pub(super) fn x_table_len(&self) -> usize {
        self.diagonals.entries.len()
    }

    /// The cell `index` of the diagonal table as `(last_hit, flag)`, for tests.
    #[cfg(test)]
    pub(super) fn x_cell(&self, index: usize) -> (i32, bool) {
        use super::x_tbn_hitbf::XDiagCell;
        match &self.diagonals.entries {
            XDiagCells::Packed(cells) => cells[index].get(),
            XDiagCells::Pair(cells) => cells[index],
        }
    }

    /// Where the diagonal table's cells are, for tests.
    #[cfg(test)]
    pub(super) fn x_table_at(&self) -> *const u8 {
        match &self.diagonals.entries {
            XDiagCells::Packed(cells) => cells.as_ptr().cast(),
            XDiagCells::Pair(cells) => cells.as_ptr().cast(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
    // The body is the per-hit loop of the two-hit word finder (the same tests, extension and HSP
    // saving of aa_ungapped.c:516-606), run on one filled offset array.
    /// The loop body of the two-hit word finder for one filled offset array.
    #[cfg(test)]
    pub(super) fn pairs(
        &mut self,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        match self.x_ungap_matrix {
            None => self.pairs_with(
                XBlosum62Walks,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
            Some(table) => {
                self.pairs_with(table, frame, chunk_offset, subject_sequence, pairs, hits)
            }
        }
    }

    /// `pairs` with the walks `walks`.
    #[inline]
    fn pairs_with<W: XUngapWalks>(
        &mut self,
        walks: W,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        // The cells are taken out for the loop over them and put back (a move of the vector).
        let mut entries =
            std::mem::replace(&mut self.diagonals.entries, XDiagCells::Pair(Vec::new()));
        match &mut entries {
            XDiagCells::Packed(cells) => self.pairs_on(
                walks,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
            XDiagCells::Pair(cells) => self.pairs_on(
                walks,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
        }
        self.diagonals.entries = entries;
    }

    /// `pairs` on the table's `cells` (either layout, `x_tbn_hitbf::XDiagCell`) with the walks
    /// `walks`. Always inlined: `x_pairs_bf_on` calls it for one pair at a time, and a call of its
    /// own set up the finder's fields for every such pair (78 instructions in a `tblastn` run of
    /// 2.26 million such calls).
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    fn pairs_on<C: super::x_tbn_hitbf::XDiagCell, W: XUngapWalks>(
        &self,
        walks: W,
        cells: &mut [C],
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        let window = self.window;
        let word_size = self.word_size;
        let word_size_i32 = word_size as i32;
        let diagonals = &self.diagonals;
        for pair in pairs {
            let q = pair.q_off as i32;
            let s = pair.s_off as i32;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-547
            let index = (q.wrapping_sub(s) & diagonals.mask) as usize;
            let cell = &mut cells[index];
            let (last_hit, flag) = cell.get();
            if flag {
                let s_diag = s.wrapping_add(diagonals.offset);
                if s_diag < last_hit {
                    continue;
                }
                cell.set(s_diag, false);
                continue;
            }
            let previous = last_hit.wrapping_sub(diagonals.offset);
            let diff = s.wrapping_sub(previous);
            if diff >= window {
                cell.set_last_hit(s.wrapping_add(diagonals.offset));
                continue;
            }
            if diff < word_size_i32 {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560-579:
            // Unsigned comparison matches NCBI's Uint4 query offset check when diff exceeds q.
            let context = self.context_offsets.partition_point(|&offset| offset <= q) - 1;
            if (q.wrapping_sub(diff) as u32) < self.context_offsets[context] as u32 {
                cell.set_last_hit(s.wrapping_add(diagonals.offset));
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-583,1089-1155:
            let s_left_off = previous.wrapping_add(word_size_i32);
            let reference = || {
                extend_two_hit(
                    self.matrix,
                    self.query_sequence,
                    subject_sequence,
                    s_left_off,
                    s as usize,
                    q as usize,
                    self.x_dropoffs[context],
                    word_size,
                )
            };
            // Only a second hit before the chunk (`s` negative, masks stored out of order with
            // `-db_soft_mask`) needs the letters around the chunk (`x_two_hit_around`). When `s`
            // is inside the chunk and `s_left_off` is negative (also reached after overlapping
            // hard masks where `x_around` is none), `extend_two_hit` and the fast walks handle
            // negative `s_left_off` directly on `subject_sequence`.
            // EXPERIMENT (LOSAT_X_TBNUNGAP): the byte-table walks (`x_tbn_ungap`); `None` = an
            // extension outside the table's range, which runs the reference.
            let result = if s < 0 {
                self.x_two_hit_around(walks, s_left_off, s, q, self.x_dropoffs[context])
            } else if self.x_ungap_fast {
                walks
                    .two_hit(
                        self.query_sequence,
                        subject_sequence,
                        s_left_off,
                        s as usize,
                        q as usize,
                        self.x_dropoffs[context],
                        word_size,
                    )
                    .unwrap_or_else(reference)
            } else {
                reference()
            };
            let Some(result) = result else {
                continue;
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590:
            if result.ungapped_data.score >= self.cutoff_scores[context] {
                hits.push(InitHsp {
                    frame,
                    chunk_offset,
                    q_seed: pair.q_off,
                    s_seed: pair.s_off,
                    q_start: result.ungapped_data.q_start,
                    s_start: result.ungapped_data.s_start,
                    length: result.ungapped_data.length,
                    score: result.ungapped_data.score,
                });
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-607
            if result.right_extend {
                cell.set(
                    result
                        .s_last_off
                        .wrapping_sub(word_size_i32 - 1)
                        .wrapping_add(diagonals.offset),
                    true,
                );
            } else {
                cell.set_last_hit(s.wrapping_add(diagonals.offset));
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-552 (the tests of
    // `pairs`; tblastn/x_tbn_hitbf.rs module doc)
    // `x_hits_bf` takes the too-far and overlap outcomes; the
    // body of `pairs` runs on each hit it stops at (a flag-1 cell or a candidate), in order.
    /// `pairs` with the too-far and overlap outcomes taken by `x_tbn_hitbf::x_hits_bf`; `pairs`
    /// itself when `x_tbn_hitbf::fits` does not hold. BLOSUM62's walks: the caller takes
    /// `x_pairs_bf_of` for a finder with another matrix's table (`x_ungap_matrix`).
    #[inline(never)]
    pub(super) fn x_pairs_bf(
        &mut self,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        // BLOSUM62's walks: never a finder with another matrix's table (`x_tbn_ungap`).
        debug_assert!(
            self.x_ungap_matrix.is_none()
                && (!self.x_ungap_fast || self.matrix == ScoringMatrix::Blosum62)
        );
        if !super::x_tbn_hitbf::fits(self.window, self.word_size) {
            self.pairs_with(
                XBlosum62Walks,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            );
            return;
        }
        // The cells are taken out for the call and put back (a move of the vector), so that the
        // layout is matched once per call.
        let mut entries =
            std::mem::replace(&mut self.diagonals.entries, XDiagCells::Pair(Vec::new()));
        match &mut entries {
            XDiagCells::Packed(cells) => self.x_pairs_bf_on(
                XBlosum62Walks,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
            XDiagCells::Pair(cells) => self.x_pairs_bf_on(
                XBlosum62Walks,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
        }
        self.diagonals.entries = entries;
    }

    /// `x_pairs_bf` with another matrix's walks, `table` (the finder's `x_ungap_matrix`). Its own
    /// text, as `x_pairs_bf`'s, so that BLOSUM62's keeps its code.
    #[inline(never)]
    pub(super) fn x_pairs_bf_of(
        &mut self,
        table: &'static x_ungapped_fast::MatrixTable,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        // The walks of the finder's own table, which `x_tbn_ungap` built for its matrix.
        debug_assert!(self.x_ungap_matrix.is_some_and(|t| std::ptr::eq(t, table)));
        if !super::x_tbn_hitbf::fits(self.window, self.word_size) {
            self.pairs_with(table, frame, chunk_offset, subject_sequence, pairs, hits);
            return;
        }
        let mut entries =
            std::mem::replace(&mut self.diagonals.entries, XDiagCells::Pair(Vec::new()));
        match &mut entries {
            XDiagCells::Packed(cells) => self.x_pairs_bf_on(
                table,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
            XDiagCells::Pair(cells) => self.x_pairs_bf_on(
                table,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
        }
        self.diagonals.entries = entries;
    }

    /// `x_pairs_bf_of` with the walks `walks` of a PSSM query (`score::XPssmWalks`, whose
    /// extension is the PSSM's in every case: the finder is made with `x_with_walks`). Its own
    /// text, as `x_pairs_bf_of`'s, so that the matrices' keep their code.
    #[inline(never)]
    pub(super) fn x_pairs_bf_walks<W: XUngapWalks>(
        &mut self,
        walks: W,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        debug_assert!(self.x_ungap_fast && self.x_ungap_matrix.is_none());
        if !super::x_tbn_hitbf::fits(self.window, self.word_size) {
            self.pairs_with(walks, frame, chunk_offset, subject_sequence, pairs, hits);
            return;
        }
        let mut entries =
            std::mem::replace(&mut self.diagonals.entries, XDiagCells::Pair(Vec::new()));
        match &mut entries {
            XDiagCells::Packed(cells) => self.x_pairs_bf_on(
                walks,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
            XDiagCells::Pair(cells) => self.x_pairs_bf_on(
                walks,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
        }
        self.diagonals.entries = entries;
    }

    /// The finder whose extensions are always its walks' (`pairs_on` takes the walks when
    /// `x_ungap_fast`): a PSSM query's (`x_pairs_bf_walks`), whatever `LOSAT_X_TBNUNGAP` and the
    /// matrix say. No NCBI counterpart: the PSSM's extension is NCBI's `use_pssm` branch, the
    /// reference of that query (`pssm_search::ungapped::pssm_extend_two_hit`).
    pub(super) fn x_with_walks(mut self) -> Self {
        self.x_ungap_fast = true;
        self.x_ungap_matrix = None;
        self
    }

    /// Another matrix's table of the byte-table walks (`x_tbn_ungap`): the caller of the word
    /// finder takes `x_pairs_bf_of` with it; `None` for BLOSUM62's walks and the reference.
    pub(super) fn x_ungap_matrix(&self) -> Option<&'static x_ungapped_fast::MatrixTable> {
        self.x_ungap_matrix
    }

    /// `x_pairs_bf` on the table's `cells` (either layout, `x_tbn_hitbf::XDiagCell`) with the
    /// walks `walks`; the loop `x_tbn_hitbf::x_hits_steps` unless `LOSAT_X_TBNHITSTEP=0` (read once
    /// per call), else `x_tbn_hitbf::x_hits_bf`.
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    fn x_pairs_bf_on<C: super::x_tbn_hitbf::XDiagCell, W: XUngapWalks>(
        &self,
        walks: W,
        cells: &mut [C],
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        if super::x_tbn_hitbf::x_hit_steps() {
            self.x_pairs_bf_with::<C, W, true>(
                walks,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            )
        } else {
            self.x_pairs_bf_with::<C, W, false>(
                walks,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            )
        }
    }

    /// `x_pairs_bf_on` with the loop `x_tbn_hitbf::x_hits::<C, STEPS>`.
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    fn x_pairs_bf_with<C: super::x_tbn_hitbf::XDiagCell, W: XUngapWalks, const STEPS: bool>(
        &self,
        walks: W,
        cells: &mut [C],
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        let word_size = self.word_size as i32;
        let mut done = 0;
        while done < pairs.len() {
            done += super::x_tbn_hitbf::x_hits::<C, STEPS>(
                cells,
                self.diagonals.offset,
                self.window,
                word_size,
                &pairs[done..],
            );
            if let Some(pair) = pairs.get(done) {
                self.pairs_on(
                    walks,
                    cells,
                    frame,
                    chunk_offset,
                    subject_sequence,
                    std::slice::from_ref(pair),
                    hits,
                );
                done += 1;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-309
    // Same order as Blast_InitHitListSortByScore: score descending, then s_start, length
    // descending, q_start. The sort is stable, as glibc's qsort is under the pinned NCBI build.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:773-803
    // The per-hit loop of `s_BlastAaWordFinder_OneHit` (`-window_size 0`; `BlastAaWordFinder`
    // takes it when the table has no `multiple_hits`, aa_ungapped.c:214-224) on one filled offset
    // array: the diagonal is `s - q` (the two-hit loop's is `q - s`), no flag is read or set, and
    // `last_hit` is set after every extension, also one below the cutoff. The table of a window of
    // 0 is never packed (`x_packs`: `x_tbn_hitbf::fits` needs a window above the word size), and
    // the unit ends with the same `exit` (`offset += subject_length + 0`). The extension is the
    // score source's (`TbnScore::one_hit`: the matrix's, BLASTP's `extend_one_hit_for`, or the
    // PSSM's `pssm_extend_one_hit`).
    /// The loop body of the one-hit word finder for one filled offset array.
    pub(super) fn x_pairs_one_hit<S: super::score::TbnScore>(
        &mut self,
        score: S,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) -> anyhow::Result<()> {
        debug_assert_eq!(self.window, 0);
        // The cells are taken out for the loop over them and put back (a move of the vector).
        let mut entries =
            std::mem::replace(&mut self.diagonals.entries, XDiagCells::Pair(Vec::new()));
        let around = self.x_around.is_some();
        let done = match &mut entries {
            XDiagCells::Packed(cells) if around => self.one_hit_on::<_, _, true>(
                score,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
            XDiagCells::Packed(cells) => self.one_hit_on::<_, _, false>(
                score,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
            XDiagCells::Pair(cells) if around => self.one_hit_on::<_, _, true>(
                score,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
            XDiagCells::Pair(cells) => self.one_hit_on::<_, _, false>(
                score,
                cells,
                frame,
                chunk_offset,
                subject_sequence,
                pairs,
                hits,
            ),
        };
        self.diagonals.entries = entries;
        done
    }

    /// `x_pairs_one_hit` on the table's `cells` (either layout, `x_tbn_hitbf::XDiagCell`).
    /// `AROUND`: the unit has the letters around its chunk (`x_set_around`), the only unit whose
    /// words can lie before its chunk; the other units' loop has no test for them.
    #[allow(clippy::too_many_arguments)]
    fn one_hit_on<
        S: super::score::TbnScore,
        C: super::x_tbn_hitbf::XDiagCell,
        const AROUND: bool,
    >(
        &self,
        score: S,
        cells: &mut [C],
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) -> anyhow::Result<()> {
        let word_size = self.word_size;
        let word_size_i32 = word_size as i32;
        let diagonals = &self.diagonals;
        for pair in pairs {
            // Uint4 offsets: the differences wrap as NCBI's unsigned arithmetic does, and the
            // Int4 `diff` is their bits.
            let q = pair.q_off as i32;
            let s = pair.s_off as i32;
            let index = (s.wrapping_sub(q) & diagonals.mask) as usize;
            let cell = &mut cells[index];
            let (last_hit, _) = cell.get();
            let diff = s.wrapping_sub(last_hit.wrapping_sub(diagonals.offset));
            if diff < 0 {
                continue;
            }
            let context = self.context_offsets.partition_point(|&offset| offset <= q) - 1;
            // The scanner gives words inside the unit or, for masks stored out of order, before
            // its chunk (the extension on the letters around it), so the extension always runs
            // (`None` is a word past the slices, which NCBI's scan does not give).
            debug_assert!(
                AROUND || s >= 0,
                "a word before a unit without letters before it"
            );
            let result = if AROUND && s < 0 {
                self.x_one_hit_around(score, q, s, self.x_dropoffs[context])
            } else {
                score.one_hit(
                    self.query_sequence,
                    subject_sequence,
                    q as usize,
                    s as usize,
                    self.x_dropoffs[context],
                    word_size,
                    self.matrix,
                )
            };
            let Some(result) = result else {
                return Err(Self::x_word_past_letters(pair.q_off, pair.s_off));
            };
            if result.ungapped_data.score >= self.cutoff_scores[context] {
                hits.push(InitHsp {
                    frame,
                    chunk_offset,
                    q_seed: pair.q_off,
                    s_seed: pair.s_off,
                    q_start: result.ungapped_data.q_start,
                    s_start: result.ungapped_data.s_start,
                    length: result.ungapped_data.length,
                    score: result.ungapped_data.score,
                });
            }
            cell.set_last_hit(result.s_last_off - (word_size_i32 - 1) + diagonals.offset);
        }
        Ok(())
    }

    // No NCBI counterpart: LOSAT's internal error for a word past the letters, which NCBI's scan
    // does not give. Out of line and by value, so that the loop keeps no copy of a pair for it.
    #[cold]
    #[inline(never)]
    fn x_word_past_letters(q_off: u32, s_off: u32) -> anyhow::Error {
        anyhow::anyhow!(
            "TBLASTN -window_size 0: the word at query offset {q_off} and subject offset {s_off} is past the letters (internal error)"
        )
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:273-313 (stable, as glibc's
    /// qsort under the pinned NCBI BLAST+).
    pub(super) fn sort(hits: &mut [InitHsp]) {
        hits.sort_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then(a.s_start.cmp(&b.s_start))
                .then(b.length.cmp(&a.length))
                .then(a.q_start.cmp(&b.q_start))
        });
    }
}

#[cfg(test)]
#[path = "search_init_tests.rs"]
mod tests;
