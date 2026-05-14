//! Optional TBLASTX stage snapshots for native/Wasm parity diagnostics.
//!
//! Set `LOSAT_DUMP_TBLASTX_STAGE=/tmp/losat-stage` to append TSV snapshots
//! under that directory. The normal runtime path is unchanged when the
//! variable is unset.

use super::chaining::UngappedHit;
use crate::common::{score_compare_hsps, Hit};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::{create_dir_all, metadata, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

struct StageDump {
    dir: PathBuf,
    lock: Mutex<()>,
}

static STAGE_DUMP: OnceLock<Option<StageDump>> = OnceLock::new();

fn stage_dump() -> Option<&'static StageDump> {
    STAGE_DUMP
        .get_or_init(|| {
            std::env::var_os("LOSAT_DUMP_TBLASTX_STAGE").map(|dir| StageDump {
                dir: PathBuf::from(dir),
                lock: Mutex::new(()),
            })
        })
        .as_ref()
}

#[inline]
pub(crate) fn enabled() -> bool {
    stage_dump().is_some()
}

fn stage_path(dump: &StageDump, stage: &str) -> PathBuf {
    let safe_stage = stage
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect::<String>();
    dump.dir.join(format!("{safe_stage}.tsv"))
}

fn write_header_if_needed(file: &mut std::fs::File, path: &PathBuf) -> std::io::Result<()> {
    if metadata(path).map(|m| m.len()).unwrap_or(0) == 0 {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:476-480
        // ```c
        // link_hsp_array =
        //    (LinkHSPStruct**) malloc(total_number_of_hsps*sizeof(LinkHSPStruct*));
        // for (index = 0; index < total_number_of_hsps; ++index) {
        //    link_hsp_array[index]->hsp = hsp_array[index];
        // ```
        writeln!(
            file,
            "stage\tq_idx\ts_idx\tctx_idx\ts_f_idx\tq_frame\ts_frame\tq_aa_start\tq_aa_end\ts_aa_start\ts_aa_end\tq_seed_off\ts_seed_off\traw_score\te_value\te_value_bits\thsp_list_order\tordering_method\tlinked_set\tstart_of_chain\tlink_id\tchain_next_link_id"
        )?;
    }
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1492-1497
// ```c
// if (Blast_HSPListReevaluateUngapped(..., hsp_list) != 0) {
//     Blast_HSPListFree(hsp_list);
//     continue;
// }
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:990-1059
// ```c
// qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//       s_RevCompareHSPsTransl);
// qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//       s_FwdCompareHSPsTransl);
// ```
// Snapshots record the internal BlastHSP/LinkHSP fields at those NCBI stage
// boundaries, before TBLASTX nucleotide-coordinate output adjustment can hide
// negative-frame or frame-relative ordering differences.
pub(crate) fn dump_ungapped_hits(stage: &str, hits: &[UngappedHit]) {
    let Some(dump) = stage_dump() else {
        return;
    };

    let _guard = dump.lock.lock().unwrap_or_else(|e| e.into_inner());
    if let Err(err) = create_dir_all(&dump.dir) {
        eprintln!(
            "[TBLASTX_STAGE_DUMP] failed to create {}: {}",
            dump.dir.display(),
            err
        );
        return;
    }

    let path = stage_path(dump, stage);
    let mut file = match OpenOptions::new().create(true).append(true).open(&path) {
        Ok(file) => file,
        Err(err) => {
            eprintln!(
                "[TBLASTX_STAGE_DUMP] failed to open {}: {}",
                path.display(),
                err
            );
            return;
        }
    };

    if let Err(err) = write_header_if_needed(&mut file, &path) {
        eprintln!(
            "[TBLASTX_STAGE_DUMP] failed to write header {}: {}",
            path.display(),
            err
        );
        return;
    }

    for h in hits {
        let chain_next = h
            .chain_next_link_id
            .map(|id| id.to_string())
            .unwrap_or_else(|| "NA".to_string());
        if let Err(err) = writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.17e}\t0x{:016x}\t{}\t{}\t{}\t{}\t{}\t{}",
            stage,
            h.q_idx,
            h.s_idx,
            h.ctx_idx,
            h.s_f_idx,
            h.q_frame,
            h.s_frame,
            h.q_aa_start,
            h.q_aa_end,
            h.s_aa_start,
            h.s_aa_end,
            h.q_seed_off,
            h.s_seed_off,
            h.raw_score,
            h.e_value,
            h.e_value.to_bits(),
            h.hsp_list_order,
            h.ordering_method,
            h.linked_set as u8,
            h.start_of_chain as u8,
            h.link_id,
            chain_next
        ) {
            eprintln!(
                "[TBLASTX_STAGE_DUMP] failed to write row {}: {}",
                path.display(),
                err
            );
            return;
        }
    }
}

#[derive(Clone)]
struct FinalSubjectGroup {
    s_idx: u32,
    best_evalue: f64,
    best_score: i32,
    hits: Vec<(Hit, UngappedHit)>,
}

#[inline]
fn evalue_comp(evalue1: f64, evalue2: f64) -> Ordering {
    const EPSILON: f64 = 1.0e-180;
    if evalue1 < EPSILON && evalue2 < EPSILON {
        Ordering::Equal
    } else if evalue1 < evalue2 {
        Ordering::Less
    } else if evalue1 > evalue2 {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3194-3205
// ```c
// Blast_HSPResultsSortByEvalue(results);
// s_Heapify((char*)hit_list->hsplist_array, ...,
//          sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1377
// ```c
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       ScoreCompareHSPs);
// ```
// Replays the same final ordering used by common.rs while dumping the paired
// internal TBLASTX HSP rows instead of formatted nucleotide coordinates.
pub(crate) fn dump_final_output_order(stage: &str, pairs: &[(Hit, UngappedHit)]) {
    if !enabled() {
        return;
    }
    if pairs.is_empty() {
        dump_ungapped_hits(stage, &[]);
        return;
    }

    let mut query_order = Vec::new();
    let mut seen_queries = HashSet::new();
    let mut query_subject_hits: HashMap<(u32, u32), Vec<(Hit, UngappedHit)>> = HashMap::new();

    for (hit, internal) in pairs.iter().cloned() {
        if seen_queries.insert(hit.q_idx) {
            query_order.push(hit.q_idx);
        }
        query_subject_hits
            .entry((hit.q_idx, hit.s_idx))
            .or_default()
            .push((hit, internal));
    }

    let mut sorted_internal = Vec::with_capacity(pairs.len());
    for q_idx in query_order {
        let subject_indices: Vec<u32> = query_subject_hits
            .keys()
            .filter(|(query_idx, _)| *query_idx == q_idx)
            .map(|(_, subject_idx)| *subject_idx)
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();

        let mut groups: Vec<FinalSubjectGroup> = subject_indices
            .into_iter()
            .filter_map(|s_idx| {
                let key = (q_idx, s_idx);
                query_subject_hits.remove(&key).map(|mut hits| {
                    hits.sort_by(|a, b| score_compare_hsps(&a.0, &b.0));
                    let best_evalue = hits
                        .iter()
                        .map(|(hit, _)| hit.e_value)
                        .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                        .unwrap_or(f64::MAX);
                    let best_score = hits.iter().map(|(hit, _)| hit.raw_score).max().unwrap_or(0);
                    FinalSubjectGroup {
                        s_idx,
                        best_evalue,
                        best_score,
                        hits,
                    }
                })
            })
            .collect();

        groups.sort_by(|a, b| {
            evalue_comp(a.best_evalue, b.best_evalue)
                .then_with(|| b.best_score.cmp(&a.best_score))
                .then_with(|| b.s_idx.cmp(&a.s_idx))
        });

        for group in groups {
            sorted_internal.extend(group.hits.into_iter().map(|(_, internal)| internal));
        }
    }

    dump_ungapped_hits(stage, &sorted_internal);
}
