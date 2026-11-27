//! EXPERIMENT (LOSAT_X_BPQBLOCK): the subjects' tabular IDs and titles of the blastp block search
//! (`x_query_blocks.rs`) made once for all blocks.
//!
//! Without these values every block's search (`run_resolved_in_pool`) makes one title `Arc` per
//! subject before its lookup table, and its output stage makes every subject's tabular ID
//! (`tabular_subject_ids`: `report::defline::tabular_subject_id`, which for a `Subject_` first word
//! generates NCBI's definition line), and both are dropped at the end of the block. These stages
//! run on the block's thread while the pool's other threads wait: in the human x yeast blastp-fast
//! job at 16 threads (workstation, `LOSAT_X_BPQBLOCK=50000`), 0.07 s for the titles, 0.16 s for the
//! IDs and part of the 0.13 s of drops at the blocks' ends, of 6.9 s (perf-mb-20261014
//! `logs/split2_report.md`).
//!
//! The block search's subjects (`x_query_blocks::XBlockSubjects::new`, made
//! once per search with the subjects' encoding) also hold an `XSubjectMeta`: the shown IDs
//! (`tabular_subject_ids`: each subject's `tabular_subject_id`, `None` where NCBI's `GenerateDefline`
//! would throw) and the title `Arc`s. A block borrows the titles (`subject_titles` in `run_resolved_in_pool`) and
//! the shown IDs, which the rows' writers (and their stops where NCBI has no result,
//! `report::tabular::row_stop_end`) read over the block's hit lists as they read
//! `tabular_subject_ids`'.
//!
//! Why the reports do not change: `tabular_subject_id(local_id, title, true)` and
//! `Arc::from(title)` are functions of the subject record only (its local ID and title bytes), and
//! the stored values are those of the same records (`x_query_blocks::x_block_subject_meta` checks
//! the slice's address and count, as `x_block_subjects` does for the encoding), so every block
//! reads the values it would make. No floating-point operation is involved. In the block search
//! no subject's ID fails: the block search is taken only when none does (`x_query_block_plan`).
//!
//! Where it acts: only the `-subject` block search, for blastp and blastp-fast (both lookup
//! tables). The search of all queries (no `LOSAT_X_BPQBLOCK`, or a search the block plan refuses)
//! attaches no block subjects to its pool, so it makes its titles and IDs as before. With `-db` the
//! search has no subject records (the database's definition lines are read per block,
//! `database::ListedTitles`), so nothing is made. `LOSAT_X_BPSUBJRUN` (the preliminary loop in runs
//! of subjects) and `LOSAT_X_BPLKPFLAT` (the lookup tables, #221) read and write none of these
//! values.
//!
//! Memory: the IDs and titles of the subject set now live for the whole block search, alongside the
//! search's peak (about 80 B per subject: about 0.5 MiB for 6,600 subjects, 8 MiB per 100k), and are
//! freed when the block search ends.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:218-227
//! NCBI makes a subject's ID from its Bioseq each time a row is written; the ID depends on the
//! record only, so one value per subject serves every block. Reuse only.

use std::sync::Arc;

use crate::blastinput::fasta_reader::FastaRecord;

/// EXPERIMENT (LOSAT_X_BPQBLOCK): the subjects' shown tabular IDs and titles, by subject index,
/// made once for all blocks of a block search.
pub(crate) struct XSubjectMeta {
    /// The IDs as `tabular_subject_ids` shows them: `tabular_subject_id(local_id, title, true)` of
    /// each subject, `None` where `GenerateDefline` throws.
    shown: Vec<crate::algorithm::common::tabular::SubjectId>,
    /// `Arc::from(record.title())` of each subject, as `run_resolved_in_pool` makes them.
    titles: Vec<Arc<[u8]>>,
}

impl XSubjectMeta {
    // No NCBI counterpart: values NCBI makes per row, made once (module comment above).
    /// The values of `records`. The block plan (`x_query_block_plan`) has already computed each ID
    /// once; this second pass (about 1-2 ms for 6,600 subjects) is kept so that the plan's hunks
    /// stay unchanged.
    pub(crate) fn new(records: &[FastaRecord]) -> Self {
        let shown = super::tabular_subject_ids(records);
        let titles = records
            .iter()
            .map(|record| Arc::from(record.title()))
            .collect();
        Self { shown, titles }
    }

    // No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
    /// The values of the block search's subject `records` when there are records (a `-db` search
    /// has none), else `None`.
    pub(crate) fn x_for_block_search(records: &[FastaRecord]) -> Option<Self> {
        (!records.is_empty()).then(|| Self::new(records))
    }

    /// The shown IDs, as `tabular_subject_ids` returns them.
    pub(crate) fn shown(&self) -> &[crate::algorithm::common::tabular::SubjectId] {
        &self.shown
    }

    /// The titles, as the reports write them.
    pub(crate) fn titles(&self) -> &[Arc<[u8]>] {
        &self.titles
    }
}

#[cfg(test)]
#[path = "x_subject_meta_tests.rs"]
mod tests;
