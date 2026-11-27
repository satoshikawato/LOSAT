//! NCBI's diagonal tables of the query batches of a protein search that holds several batches in
//! one search (BLASTP's search of all queries, or a block of whole batches of LOSAT_X_BPQBLOCK),
//! for the `-db_hard_mask` sequences whose hits can meet entries that the scans of their own
//! chunks did not write (series NP, M4b-3b): [`QueryBatchTables`].
//!
//! NCBI searches each query batch on its own (`blastp_app.cpp:216, 257`), on a table of its own
//! whose diagonals and query offsets count from the batch's start. A search of several batches at
//! once on one table gives each batch NCBI's result only where every entry within reach of a hit
//! was written by a hit of the same scan, forward over the sequence: on each diagonal the hits of a
//! batch then follow those of the batches before it, and an entry of another batch reads as one
//! out of reach in both (`blastp::x_query_blocks`, the header). An entry that no sequence wrote
//! (a new table at NCBI's one-thread offset, `common::hard_chunk_walk::TableStart::New`) or one of
//! an earlier chunk of the sequence after a walk back reads otherwise: the batch's own hits on
//! that diagonal, NCBI's ring of the batch's size and NCBI's context test of the batch's offsets
//! decide, which one table of all batches cannot hold. Such a sequence is searched on these
//! tables, every hit on its own batch's ring (`common::hard_chunk_walk::DiagonalWalk::apart`).

use super::diag_entry::{DiagStruct, DiagTable};
use super::hard_chunk_walk::{HardTable, NcbiStart};

/// NCBI's ring of one query batch: its entries (made at the first hit of a sequence searched on
/// the tables, `mark`) and its mask.
#[derive(Debug, Clone)]
struct BatchRing {
    entries: Vec<DiagStruct>,
    mask: i32,
    mark: u64,
}

/// NCBI's diagonal table of each query batch of a search, all at one diagonal offset (NCBI's
/// offset after a sequence depends only on the sequences: the same in every batch's search),
/// each ring made as the table of the sequence's start where a hit of the sequence first reads it
/// ([`HardTable::ncbi_table`]).
#[derive(Debug, Clone)]
pub(crate) struct QueryBatchTables {
    /// Each batch's first offset in the search's concatenated query and its concatenated length
    /// (the layout of [`QueryBatchTables::new`]).
    layout: Vec<(i32, i32)>,
    /// The search's window (the rings' sizes hold it).
    window: i32,
    rings: Vec<BatchRing>,
    /// The diagonal offset of every ring.
    pub(crate) offset: i32,
    /// The entry of a ring that no hit of the sequence has read yet.
    fill: DiagStruct,
    /// The mark of the sequence's tables (a ring whose mark is older holds another's entries).
    mark: u64,
}

impl QueryBatchTables {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1004
    /// The tables of the query batches `layout` (each batch's first offset in the concatenated
    /// query and its concatenated length, in order; one batch of all queries where it is empty)
    /// in a search with the window `window`: each ring the power of two just at least the
    /// batch's length plus the window, as NCBI makes the table of a batch's search.
    pub(crate) fn new(layout: &[(i32, i32)], window: i32) -> Self {
        let layout = if layout.is_empty() {
            vec![(0, 0)]
        } else {
            layout.to_vec()
        };
        let rings = layout
            .iter()
            .map(|&(_, length)| {
                let mut size = 1i64;
                while size < i64::from(length) + i64::from(window) {
                    size <<= 1;
                }
                BatchRing {
                    entries: Vec::new(),
                    mask: i32::try_from(size - 1).unwrap_or(i32::MAX),
                    mark: 0,
                }
            })
            .collect();
        Self {
            layout,
            window,
            rings,
            offset: window,
            fill: DiagStruct::default(),
            mark: 1,
        }
    }

    /// Whether these are the tables of the batches `layout` with the window `window` (a worker's
    /// tables are kept for the search that made them).
    pub(crate) fn is_of(&self, layout: &[(i32, i32)], window: i32) -> bool {
        self.window == window
            && if layout.is_empty() {
                self.layout == [(0, 0)]
            } else {
                self.layout == layout
            }
    }

    /// The ring of the batch of query offset `q` (made as the table of the sequence's start where
    /// the sequence has not read it), and the batch's first offset.
    #[inline]
    fn ring_of(&mut self, q: i32) -> (&mut BatchRing, i32) {
        let index = self
            .layout
            .partition_point(|&(start, _)| start <= q)
            .saturating_sub(1);
        let ring = &mut self.rings[index];
        if ring.mark != self.mark {
            let size = ring.mask as usize + 1;
            if ring.entries.len() == size {
                ring.entries.fill(self.fill);
            } else {
                ring.entries = vec![self.fill; size];
            }
            ring.mark = self.mark;
        }
        (ring, self.layout[index].0)
    }

    /// The table of a new sequence, every entry `fill`, at the offset `offset`.
    fn restart(&mut self, fill: DiagStruct, offset: i32) {
        self.fill = fill;
        self.offset = offset;
        self.mark += 1;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
    /// `Blast_ExtendWordExit` after a chunk of `length` letters (NCBI's `Int4`), on every ring.
    pub(crate) fn finish_chunk(&mut self, length: i32, window: i32) {
        if self.offset >= i32::MAX / 4 {
            self.clear_table(window);
        } else {
            self.offset += length + window;
        }
    }

    /// The rings made so far (for the tests: a batch without hits makes none).
    #[cfg(test)]
    pub(crate) fn made(&self) -> Vec<usize> {
        self.rings
            .iter()
            .map(|ring| {
                if ring.mark == self.mark {
                    ring.entries.len()
                } else {
                    0
                }
            })
            .collect()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-518, 776-779
// The diagonal of a hit in the offsets of its batch's search.
impl DiagTable for QueryBatchTables {
    type Cell = DiagStruct;

    fn two_hit_entry(&mut self, q: i32, s: i32) -> (&mut DiagStruct, i32) {
        let (ring, start) = self.ring_of(q);
        let index = (q.wrapping_sub(start).wrapping_sub(s) & ring.mask) as usize;
        (&mut ring.entries[index], start)
    }

    fn one_hit_entry(&mut self, q: i32, s: i32) -> &mut DiagStruct {
        let (ring, start) = self.ring_of(q);
        let index = (s.wrapping_sub(q.wrapping_sub(start)) & ring.mask) as usize;
        &mut ring.entries[index]
    }
}

// The search of a sequence on the tables starts them as NCBI's table at the sequence's start
// (`common::hard_chunk_walk::DiagonalWalk::apart`); a clear inside it clears every ring.
impl HardTable for QueryBatchTables {
    type Error = std::convert::Infallible;

    fn ncbi_table(&mut self, start: NcbiStart, window: i32) {
        let fill = if start.cleared {
            DiagStruct::clear(window)
        } else {
            DiagStruct::default()
        };
        self.restart(fill, start.offset);
    }

    fn clear_table(&mut self, window: i32) {
        self.restart(DiagStruct::clear(window), window);
    }

    fn move_offset(&mut self, by: i32) -> Result<(), Self::Error> {
        self.offset += by;
        Ok(())
    }
}

#[cfg(test)]
#[path = "query_batch_tables_tests.rs"]
mod tests;
