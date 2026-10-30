//! BLASTP's Seq-annot of Seq-aligns (outfmt 8, 12): the report records of its final hit lists
//! (`report::record`), one query at a time, written by `report::serial::write_query`.
//!
//! NCBI makes the Seq-aligns of each query from its final hit list
//! (`BlastHitList2SeqAlign_OMF`, see `algorithm::blastn::seqalign`); BLASTP's reports are made
//! from the lists after `hit_list_to_seqalign_order`, which puts them in that order, so the
//! documents walk them as the tabular rows do (`write_blastp_tabular_hit_lists`).
//!
//! What NCBI's `BlastHSP` holds, from LOSAT's `BlastpHsp` (after `shift_hsp_to_records`,
//! which moved the coordinates to the records):
//! - `query.offset/end`: `q_start - 1` and `q_end` less the `-query_loc` start;
//!   `subject.offset/end`: `s_start - 1` and `s_end` less the `-subject_loc` start
//!   (`RemapToQueryLoc` / `s_RemapToSubjectLoc` add them back to the rows);
//! - frames 0 (proteins: strand `unknown`); `gap_info`: `gap_info`;
//! - `score`, `evalue`, `bit_score`, `comp_adjustment_method`: `raw_score`, `e_value`,
//!   `bit_score`, `comp_adjustment_method`; `num` 0 (no linking);
//! - `num_ident`, `num_positives`: the engine's counts (`kappa.rs` `recompute_num_identities`,
//!   NCBI's `Blast_HSPGetNumIdentitiesAndPositives` on the search's letters, where O is read
//!   as X), which `recount_shown_identities` replaces by the display counts for the pairs with O
//!   and the subjects with a gap character; `EngineCounts` keeps the engine's values of those
//!   HSPs.
//! - the query length of the construction: the searched length (`query_lengths`, with
//!   `-query_loc` the interval's), not `BlastpHsp::query_length`, which holds the range as
//!   given.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:518-531
//! ```c
//!             query = query_blk->sequence +
//!                 query_info->contexts[hsp->context].query_offset;
//!             query_nomask = query_blk->sequence_nomask +
//!                 query_info->contexts[hsp->context].query_offset;
//!         }
//!
//!         /* Translate subject if needed. */
//!         if (program_number == eBlastTypeTblastn) {
//!             ...
//!         }
//!         else
//!             status = Blast_HSPGetNumIdentitiesAndPositives(query_nomask, subject, hsp, scoring_options, 0, sbp);
//! ```

use std::io::Write;

use anyhow::Result;

use super::hsp::BlastpHitList;
use crate::blastinput::app::Program;
use crate::blastinput::seq_range::Placements;
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryNames, RunRecord, SubjectNames, SubjectRecord};
use crate::report::serial::{write_hit_list_seqaligns, SerialFormat};

/// Where an HSP lies in the final hit lists: (query, HSP list, HSP).
type HspKey = (u32, u32, u32);

/// The engine's `(num_ident, num_positives)` of the HSPs whose counts the display recount
/// replaced (`recount_shown_identities`), kept only when a Seq-annot format is
/// requested (`ReportNeeds::seqalign`). The keys are added in the order of the hit lists, so
/// they are sorted.
#[derive(Default)]
pub(crate) struct EngineCounts {
    saved: Vec<(HspKey, (usize, usize))>,
}

impl EngineCounts {
    /// Keeps the engine's counts of the HSP at `key` before they are replaced.
    pub(crate) fn save(&mut self, key: HspKey, counts: (usize, usize)) {
        debug_assert!(self.saved.last().is_none_or(|&(last, _)| last < key));
        self.saved.push((key, counts));
    }

    /// The engine's counts of the HSP at `key`, if the recount replaced them.
    fn get(&self, key: HspKey) -> Option<(usize, usize)> {
        self.saved
            .binary_search_by(|(saved, _)| saved.cmp(&key))
            .ok()
            .map(|index| self.saved[index].1)
    }
}

/// Writes the documents of every query of `hit_lists` (indexed by query) to `writer`, each
/// query's warnings before it.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_blastp_seqalign(
    format: SerialFormat,
    hit_lists: &[Option<BlastpHitList>],
    queries: &QueryNames<'_>,
    subjects: SubjectNames<'_>,
    subject_placements: &Placements,
    engine_counts: &EngineCounts,
    writer: &mut dyn Write,
    warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Blastp,
        db_title: subjects.db_title(),
        gapped: true,
    };
    write_hit_list_seqaligns(
        format,
        &run,
        hit_lists,
        queries,
        writer,
        warnings,
        |q_idx, hit_list, subject_records, hsp_records| {
            query_records(
                hit_list,
                q_idx,
                &subjects,
                subject_placements,
                engine_counts,
                queries.placements.offset(q_idx) as i32,
                subject_records,
                hsp_records,
            );
        },
    )
}

/// The subject and HSP records of query `q_idx`'s hit list (what `write_query` and the archive
/// write), replacing the contents of `subject_records` and `hsp_records`. `q_shift` is the start
/// of the query's interval in its record (the HSPs' coordinates are the record's, less the shift).
#[allow(clippy::too_many_arguments)]
pub(crate) fn query_records<'a>(
    hit_list: &'a BlastpHitList,
    q_idx: usize,
    subjects: &SubjectNames<'a>,
    subject_placements: &Placements,
    engine_counts: &EngineCounts,
    q_shift: i32,
    subject_records: &mut Vec<SubjectRecord<'a>>,
    hsp_records: &mut Vec<HspRecord<'a>>,
) {
    subject_records.clear();
    hsp_records.clear();
    for (list_index, hsp_list) in hit_list
        .hsplist_array
        .iter()
        .take(hit_list.hsplist_count)
        .enumerate()
    {
        if hsp_list.hsps.is_empty() {
            continue;
        }
        let s_idx = hsp_list.hsps[0].s_idx as usize;
        let subject = subject_records.len() as u32;
        // A protein subject row (strand unknown) does not read its length.
        subject_records.push(subjects.record(s_idx, subject_placements.offset(s_idx), 0));
        let s_shift = subject_placements.offset(s_idx) as i32;
        for (hsp_index, hsp) in hsp_list.hsps.iter().enumerate() {
            debug_assert_eq!(hsp.s_idx as usize, s_idx, "one subject per HSP list");
            let (num_ident, num_positives) = engine_counts
                .get((q_idx as u32, list_index as u32, hsp_index as u32))
                .unwrap_or((hsp.num_ident, hsp.num_positives));
            hsp_records.push(HspRecord {
                subject,
                q_offset: hsp.q_start as i32 - 1 - q_shift,
                q_end: hsp.q_end as i32 - q_shift,
                s_offset: hsp.s_start as i32 - 1 - s_shift,
                s_end: hsp.s_end as i32 - s_shift,
                q_frame: 0,
                s_frame: 0,
                score: hsp.raw_score,
                num: 0,
                evalue: hsp.e_value,
                bit_score: hsp.bit_score,
                num_ident: num_ident as i32,
                num_positives: num_positives as i32,
                comp_adjust: hsp.comp_adjustment_method,
                ops: hsp.gap_info.as_deref().unwrap_or(&[]),
            });
        }
    }
}
