//! EXPERIMENT (LOSAT_X_BPSUBJMETA / LOSAT_X_BPSUBJMETASHADOW): the subjects' tabular IDs and titles
//! of the blastp block search (`LOSAT_X_BPQBLOCK`, `x_query_blocks.rs`) made once for all blocks.
//!
//! Without the switch every block's search (`run_resolved_in_pool`) makes one title `Arc` per
//! subject before its lookup table, and its output stage makes every subject's tabular ID
//! (`tabular_subject_ids`: `report::defline::tabular_subject_id`, which for a `Subject_` first word
//! generates NCBI's definition line), and both are dropped at the end of the block. These stages
//! run on the block's thread while the pool's other threads wait: in the human x yeast blastp-fast
//! job at 16 threads (workstation, `LOSAT_X_BPQBLOCK=50000`), 0.07 s for the titles, 0.16 s for the
//! IDs and part of the 0.13 s of drops at the blocks' ends, of 6.9 s (perf-mb-20261014
//! `logs/split2_report.md`).
//!
//! With the switch on, the block search's subjects (`x_query_blocks::XBlockSubjects::new`, made
//! once per search with the subjects' encoding) also hold an `XSubjectMeta`: each subject's
//! `tabular_subject_id` result (`None` where NCBI's `GenerateDefline` would throw; not an error
//! until a row needs it), the shown IDs (empty for `None`, as `tabular_subject_ids` shows them) and
//! the title `Arc`s. A block borrows the titles (`subject_titles` in `run_resolved_in_pool`) and
//! the shown IDs (`x_tabular_subject_ids_cached`, which runs the check loop of
//! `tabular_subject_ids` over the block's hit lists with the stored results).
//!
//! Why the reports do not change: `tabular_subject_id(local_id, title, true)` and
//! `Arc::from(title)` are functions of the subject record only (its local ID and title bytes), and
//! the stored values are those of the same records (`x_query_blocks::x_block_subject_meta` checks
//! the slice's address and count, as `x_block_subjects` does for the encoding), so every block
//! reads the values it would make. The error for a subject with a row and no decodable ID is raised
//! by the same loop (`check_tabular_subject_ids`) over the same block's hit lists with the same
//! per-subject results, so the same block fails, before any of its rows, naming the same record. No
//! floating-point operation is involved. In the block search the error does not arise at all: the
//! block search is taken only when no subject's ID fails (`x_query_block_plan`); the loop is kept
//! as the reference path's.
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
//! `LOSAT_X_BPSUBJMETA` is a switch of the strict set (on by default since 2026-10-17): on unless it is `0`.
//!
//! `LOSAT_X_BPSUBJMETASHADOW=1` (shadow mode, implies the switch): in each block whose output stage
//! takes the stored IDs, the IDs and titles are made again as without the switch and compared with
//! the stored ones element for element: each subject's `tabular_subject_id` result, its title, and
//! the result of `tabular_subject_ids` (the shown IDs, or the error's text). The first difference
//! stops the run (panic) after printing the line that is otherwise printed at exit: `[X_SHADOW]
//! BPSUBJMETA blocks= subjects_compared= identical` (`blocks`: output stages that took the stored
//! values; `subjects_compared`: the subjects compared in them). The output stage of a block's
//! tabular report (outfmt 6, or 7 with `LOSAT_X_BPQBLOCKFMT`) makes the IDs and is counted; an
//! outfmt 0 report (`LOSAT_X_BPQBLOCKFMT`) makes none.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:218-227
//! ```c++
//!     ITERATE(vector< CConstRef<CSeq_id> >, itr, original_seqids) {
//!         CRef<CSeq_id> next_seqid(new CSeq_id());
//!         string id_token = NcbiEmptyString;
//!
//!         if (((*itr)->IsGeneral() &&
//!             (*itr)->AsFastaString().find("gnl|BL_ORD_ID")
//!             != string::npos) ||
//! 		(*itr)->AsFastaString().find("lcl|Subject_") != string::npos) {
//!             vector<string> title_tokens;
//!             string defline = sequence::CDeflineGenerator().GenerateDefline(bh);
//! ```
//! NCBI makes a subject's ID from its Bioseq each time a row is written; the ID depends on the
//! record only, so one value per subject serves every block. Reuse only.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use crate::blastinput::fasta_reader::FastaRecord;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XSubjMeta {
    Off,
    On,
    Shadow,
}

// No NCBI counterpart: reads the LOSAT_X_BPSUBJMETA / LOSAT_X_BPSUBJMETASHADOW switches once; it
// does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPSUBJMETA): a strict-set switch, on unless it is `0`;
/// `LOSAT_X_BPSUBJMETASHADOW` (any value) selects the shadow
/// mode.
pub(crate) fn x_bp_subj_meta() -> XSubjMeta {
    static MODE: OnceLock<XSubjMeta> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BPSUBJMETASHADOW").is_some() {
            XSubjMeta::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_BPSUBJMETA") {
            XSubjMeta::On
        } else {
            XSubjMeta::Off
        }
    })
}

/// EXPERIMENT (LOSAT_X_BPSUBJMETA): the subjects' tabular ID results, shown IDs and titles, by
/// subject index, made once for all blocks of a block search.
pub(crate) struct XSubjectMeta {
    /// `tabular_subject_id(local_id, title, true)` of each subject
    /// (`super::tabular_subject_id_results`).
    ids: Vec<Option<Arc<[u8]>>>,
    /// The IDs as `tabular_subject_ids` shows them: the empty ID for `None`.
    shown: Vec<Arc<[u8]>>,
    /// `Arc::from(record.title())` of each subject, as `run_resolved_in_pool` makes them.
    titles: Vec<Arc<[u8]>>,
    /// The shadow mode (`LOSAT_X_BPSUBJMETASHADOW`) was on when they were made.
    shadow: bool,
}

impl XSubjectMeta {
    // No NCBI counterpart: values NCBI makes per row, made once (module comment above).
    /// The values of `records`, with the shadow mode's flag `shadow`. The block plan
    /// (`x_query_block_plan`) has already computed each ID once; this second pass (about 1-2 ms for
    /// 6,600 subjects) is kept so that the plan's hunks stay unchanged.
    pub(crate) fn new(records: &[FastaRecord], shadow: bool) -> Self {
        let ids = super::tabular_subject_id_results(records);
        let shown = ids
            .iter()
            .map(|id| id.clone().unwrap_or_else(|| Arc::from(&b""[..])))
            .collect();
        let titles = records
            .iter()
            .map(|record| Arc::from(record.title()))
            .collect();
        Self {
            ids,
            shown,
            titles,
            shadow,
        }
    }

    // No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
    /// The values of the block search's subject `records` when the switch is on and there are
    /// records (a `-db` search has none), else `None`.
    pub(crate) fn x_for_block_search(records: &[FastaRecord]) -> Option<Self> {
        let mode = x_bp_subj_meta();
        if mode != XSubjMeta::Off {
            crate::utils::xreach::x_reach!("LOSAT_X_BPSUBJMETA", !records.is_empty());
        }
        (mode != XSubjMeta::Off && !records.is_empty())
            .then(|| Self::new(records, mode == XSubjMeta::Shadow))
    }

    /// Each subject's `tabular_subject_id` result (`None`: the ID cannot be written).
    pub(crate) fn ids(&self) -> &[Option<Arc<[u8]>>] {
        &self.ids
    }

    /// The shown IDs, as `tabular_subject_ids` returns them when no reported subject fails.
    pub(crate) fn shown(&self) -> &[Arc<[u8]>] {
        &self.shown
    }

    /// The titles, as the reports write them.
    pub(crate) fn titles(&self) -> &[Arc<[u8]>] {
        &self.titles
    }

    /// Whether the shadow mode compares the values of each block.
    pub(crate) fn shadow(&self) -> bool {
        self.shadow
    }

    // No NCBI counterpart: a check of the shadow mode; it changes no value NCBI computes.
    /// EXPERIMENT (LOSAT_X_BPSUBJMETASHADOW): the stored values against those made again for a
    /// block of `records` as without the switch: `reference_ids` (`tabular_subject_id_results`),
    /// the titles, and `reference` (`tabular_subject_ids` of the block's hit lists) against
    /// `cached` (the result the block takes). Panics at the first difference.
    pub(crate) fn x_shadow_compare(
        &self,
        records: &[FastaRecord],
        reference_ids: &[Option<Arc<[u8]>>],
        reference: &anyhow::Result<Vec<Arc<[u8]>>>,
        cached: &anyhow::Result<&[Arc<[u8]>]>,
    ) {
        if self.ids.len() != records.len()
            || reference_ids.len() != records.len()
            || self.titles.len() != records.len()
            || self.shown.len() != records.len()
        {
            x_shadow_difference(format!(
                "{} records, {} reference IDs, {} stored IDs, {} shown IDs, {} titles",
                records.len(),
                reference_ids.len(),
                self.ids.len(),
                self.shown.len(),
                self.titles.len()
            ));
        }
        for (index, record) in records.iter().enumerate() {
            if self.ids[index].as_deref() != reference_ids[index].as_deref() {
                x_shadow_difference(format!(
                    "subject {index}: stored ID {:?}, reference {:?}",
                    self.ids[index].as_deref().map(String::from_utf8_lossy),
                    reference_ids[index].as_deref().map(String::from_utf8_lossy)
                ));
            }
            if *self.titles[index] != *record.title() {
                x_shadow_difference(format!(
                    "subject {index}: stored title {:?}, reference {:?}",
                    String::from_utf8_lossy(&self.titles[index]),
                    String::from_utf8_lossy(record.title())
                ));
            }
        }
        match (reference, cached) {
            (Ok(reference), Ok(cached)) => {
                if reference.len() != cached.len() {
                    x_shadow_difference(format!(
                        "{} shown IDs, reference {}",
                        cached.len(),
                        reference.len()
                    ));
                }
                for (index, (reference, cached)) in reference.iter().zip(cached.iter()).enumerate()
                {
                    if **reference != **cached {
                        x_shadow_difference(format!(
                            "subject {index}: shown ID {:?}, reference {:?}",
                            String::from_utf8_lossy(cached),
                            String::from_utf8_lossy(reference)
                        ));
                    }
                }
            }
            (Err(reference), Err(cached)) => {
                let (reference, cached) = (format!("{reference:#}"), format!("{cached:#}"));
                if reference != cached {
                    x_shadow_difference(format!("error {cached:?}, reference {reference:?}"));
                }
            }
            (reference, cached) => x_shadow_difference(format!(
                "result {:?}, reference {:?}",
                cached
                    .as_ref()
                    .map(|ids| ids.len())
                    .map_err(|e| e.to_string()),
                reference.as_ref().map(Vec::len).map_err(|e| e.to_string())
            )),
        }
        X_BLOCKS.fetch_add(1, Ordering::Relaxed);
        X_SUBJECTS_COMPARED.fetch_add(records.len() as u64, Ordering::Relaxed);
    }
}

static X_BLOCKS: AtomicU64 = AtomicU64::new(0);
static X_SUBJECTS_COMPARED: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: prints counters and changes no value NCBI computes.
fn x_shadow_print_line(identical: bool) {
    eprintln!(
        "[X_SHADOW] BPSUBJMETA blocks={} subjects_compared={}{}",
        X_BLOCKS.load(Ordering::Relaxed),
        X_SUBJECTS_COMPARED.load(Ordering::Relaxed),
        if identical {
            " identical"
        } else {
            " DIFFERENT"
        },
    );
}

// No NCBI counterpart: prints counters and changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPSUBJMETASHADOW): the line `main` prints at exit.
pub(crate) fn x_bp_subj_meta_print_summary() {
    if x_bp_subj_meta() == XSubjMeta::Shadow {
        x_shadow_print_line(true);
    }
}

// No NCBI counterpart: the shadow mode's check stops the run on a difference.
fn x_shadow_difference(message: String) -> ! {
    // The panic ends the run before `main` prints the line at exit: it is printed here.
    x_shadow_print_line(false);
    panic!("[X_SHADOW] BPSUBJMETA {message}");
}

#[cfg(test)]
mod tests {
    use super::super::{
        check_tabular_subject_ids, tabular_subject_id_results, tabular_subject_ids,
        x_tabular_subject_ids_cached,
    };
    use super::*;
    use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHsp, BlastpHspList};
    use crate::common::Hit;

    /// Subject records whose IDs take every branch of `tabular_subject_id`: a first word, the local
    /// ID of an empty title, `Subject_` first words that `GenerateDefline` keeps, rewrites or
    /// cannot decode (`\x81`, no encoding NCBI guesses: `None`), and an empty record.
    fn records() -> Vec<FastaRecord> {
        let titles: [&[u8]; 9] = [
            b"sp|P1|A first protein",
            b"",
            b"Subject_1 desc",
            b"Subject_\x81 desc",
            b"   ",
            b"Subject_\xe9 desc",
            b"lcl|Subject_7 foo",
            b"Subject_\x81",
            b"x",
        ];
        titles
            .iter()
            .enumerate()
            .map(|(index, title)| {
                FastaRecord::new(format!("Subject_{}", index + 1), title, b"MKVLAAGIVA")
            })
            .collect()
    }

    /// A hit list of query 0 with one HSP list per subject index of `subjects`.
    fn hit_list(subjects: &[u32]) -> Option<BlastpHitList> {
        let mut list = BlastpHitList::new(subjects.len().max(1));
        for &s_idx in subjects {
            let hit = Hit {
                identity: 100.0,
                length: 10,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 10,
                s_start: 1,
                s_end: 10,
                e_value: 1.0,
                bit_score: 0.0,
                num_ident: 10,
                query_frame: 0,
                query_length: 10,
                q_idx: 0,
                s_idx,
                raw_score: 50,
                sort_query_offset: 0,
                sort_query_end: 0,
                sort_subject_offset: 0,
                sort_subject_end: 0,
                has_sort_offsets: false,
                gap_info: None,
                num_positives: 10,
            };
            list.hsplist_array.push(BlastpHspList {
                oid: s_idx,
                query_index: 0,
                hsps: vec![BlastpHsp::from_hit(hit)],
                best_evalue: 1.0,
            });
        }
        list.hsplist_count = list.hsplist_array.len();
        Some(list)
    }

    // The stored values are those that every block makes: each subject's ID result (including the
    // records whose ID fails), the shown IDs of `tabular_subject_ids` and the titles.
    #[test]
    fn stored_values_are_the_per_block_values() {
        let records = records();
        let meta = XSubjectMeta::new(&records, false);
        let reference_ids = tabular_subject_id_results(&records);
        assert_eq!(meta.ids().len(), records.len());
        for (index, record) in records.iter().enumerate() {
            let direct = crate::report::defline::tabular_subject_id(
                record.local_id.as_bytes(),
                &record.title,
                true,
            )
            .ok();
            assert_eq!(meta.ids()[index].as_deref(), direct.as_deref(), "{index}");
            assert_eq!(
                meta.ids()[index].as_deref(),
                reference_ids[index].as_deref(),
                "{index}"
            );
            assert_eq!(&*meta.titles()[index], record.title(), "{index}");
        }
        // Two of the records cannot be written (`Subject_\x81`), the others can.
        let failing: Vec<usize> = (0..records.len())
            .filter(|&index| meta.ids()[index].is_none())
            .collect();
        assert_eq!(failing, vec![3, 7]);
        // No reported subject: `tabular_subject_ids` shows every ID, the empty one where it fails.
        let reference = tabular_subject_ids(&[], &records).expect("no subject is reported");
        let shown: Vec<&[u8]> = meta.shown().iter().map(|id| &id[..]).collect();
        let reference: Vec<&[u8]> = reference.iter().map(|id| &id[..]).collect();
        assert_eq!(shown, reference);
        assert_eq!(shown[3], b"");
        assert_eq!(shown[7], b"");
        assert_eq!(shown[0], b"sp|P1|A");
    }

    // The check over a block's hit lists gives the result of `tabular_subject_ids`: the shown IDs
    // when no reported subject fails, else the same error naming the same record (the first in the
    // hit lists' order), whatever other subjects fail without rows.
    #[test]
    fn the_check_gives_the_per_block_result_and_error() {
        let records = records();
        let meta = XSubjectMeta::new(&records, true);
        let cases: Vec<Vec<Option<BlastpHitList>>> = vec![
            vec![],
            vec![None, hit_list(&[0, 2, 5])],
            vec![hit_list(&[8]), None, hit_list(&[1, 6])],
            vec![hit_list(&[0, 3])],
            vec![hit_list(&[2]), hit_list(&[7, 3])],
            vec![hit_list(&[3]), hit_list(&[7])],
        ];
        let mut errors = 0;
        for (case, hit_lists) in cases.iter().enumerate() {
            let reference = tabular_subject_ids(hit_lists, &records);
            let cached = x_tabular_subject_ids_cached(hit_lists, &records, &meta);
            match (&reference, &cached) {
                (Ok(reference), Ok(cached)) => {
                    assert_eq!(reference.len(), cached.len(), "{case}");
                    assert!(
                        reference.iter().zip(cached.iter()).all(|(a, b)| a == b),
                        "{case}"
                    );
                }
                (Err(reference), Err(cached)) => {
                    errors += 1;
                    assert_eq!(format!("{reference:#}"), format!("{cached:#}"), "{case}");
                }
                _ => panic!("case {case}: {reference:?} against {cached:?}"),
            }
            // The stored IDs give the same check result as the IDs made again.
            let again = check_tabular_subject_ids(
                hit_lists,
                &records,
                &tabular_subject_id_results(&records),
            );
            assert_eq!(again.is_ok(), cached.is_ok(), "{case}");
        }
        assert_eq!(errors, 3);
        let error = x_tabular_subject_ids_cached(&cases[4], &records, &meta).unwrap_err();
        assert!(
            format!("{error:#}").starts_with("subject record 8 (Subject_"),
            "{error:#}"
        );
    }

    // The shadow's comparison stops at a stored value that differs (here: the titles and IDs of
    // other records), after printing its line.
    #[test]
    fn the_shadow_stops_at_a_difference() {
        let records = records();
        let mut others = records.clone();
        others.swap(0, 2);
        let meta = XSubjectMeta::new(&others, true);
        let reference_ids = tabular_subject_id_results(&records);
        let reference = tabular_subject_ids(&[], &records);
        let cached: anyhow::Result<&[Arc<[u8]>]> = Ok(meta.shown());
        let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            meta.x_shadow_compare(&records, &reference_ids, &reference, &cached)
        }))
        .expect_err("a difference");
        let message = panic.downcast_ref::<String>().cloned().unwrap_or_default();
        assert!(
            message.starts_with("[X_SHADOW] BPSUBJMETA subject 0: stored ID"),
            "{message}"
        );
        // The same records compare equal.
        let meta = XSubjectMeta::new(&records, true);
        let cached: anyhow::Result<&[Arc<[u8]>]> = Ok(meta.shown());
        meta.x_shadow_compare(&records, &reference_ids, &reference, &cached);
    }
}
