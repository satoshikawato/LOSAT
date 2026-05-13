//! BLASTN HSP tracing infrastructure for parity diagnostics.
//!
//! Enable with:
//! - `LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"`
//! - `LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>`
//! - `LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>`
//! - `LOSAT_TRACE_BLASTN_STAGE=<seed|ungapped|prelim|traceback|purge|hitlist|all>`
//!
//! The target coordinates use final outfmt 6/7 coordinates, matching NCBI
//! tabular output. Earlier pipeline stages are matched by overlap with the
//! target range so a selected final HSP can be followed back through seeds,
//! ungapped extension, preliminary gapped extension, traceback, and pruning.

use std::sync::OnceLock;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1109-1132
// ```c
// void Blast_HSPGetAdjustedOffsets(EBlastProgramType program, BlastHSP* hsp,
//                                  Int4 query_length, Int4* q_start,
//                                  Int4* q_end, Int4* s_start, Int4* s_end)
// {
//    if (hsp->query.frame != hsp->subject.frame) {
//       *q_end = query_length - hsp->query.offset;
//       *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
//       *s_end = hsp->subject.offset + 1;
//       *s_start = hsp->subject.end;
//    } else {
//       *q_start = hsp->query.offset + 1;
//       *q_end = hsp->query.end;
//       *s_start = hsp->subject.offset + 1;
//       *s_end = hsp->subject.end;
//    }
// }
// ```
#[derive(Clone, Copy, Debug)]
pub struct TraceBlastnTarget {
    pub q_start: i64,
    pub q_end: i64,
    pub s_start: i64,
    pub s_end: i64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum StageFilter {
    All,
    Seed,
    Ungapped,
    Prelim,
    Traceback,
    Purge,
    Hitlist,
}

#[derive(Clone, Debug)]
struct TraceBlastnConfig {
    target: Option<TraceBlastnTarget>,
    context_idx: Option<u32>,
    subject: Option<String>,
    stage: StageFilter,
    enabled: bool,
}

static TRACE_BLASTN_CONFIG: OnceLock<TraceBlastnConfig> = OnceLock::new();

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., FALSE);
// ...
// Blast_HSPReevaluateWithAmbiguitiesGapped(...);
// ...
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE);
// Blast_HSPListSortByScore(hsp_list);
// Blast_IntervalTreeReset(tree);
// ```
fn config() -> &'static TraceBlastnConfig {
    TRACE_BLASTN_CONFIG.get_or_init(|| {
        let target = std::env::var("LOSAT_TRACE_BLASTN_HSP")
            .ok()
            .and_then(|raw| parse_target(raw.as_str()));
        let context_idx = std::env::var("LOSAT_TRACE_BLASTN_CONTEXT")
            .ok()
            .and_then(|raw| raw.trim().parse::<u32>().ok());
        let subject = std::env::var("LOSAT_TRACE_BLASTN_SUBJECT")
            .ok()
            .map(|raw| raw.trim().to_owned())
            .filter(|raw| !raw.is_empty());
        let stage = std::env::var("LOSAT_TRACE_BLASTN_STAGE")
            .ok()
            .as_deref()
            .and_then(parse_stage)
            .unwrap_or(StageFilter::All);
        let enabled = target.is_some()
            || context_idx.is_some()
            || subject.is_some()
            || std::env::var("LOSAT_TRACE_BLASTN_STAGE").is_ok();

        TraceBlastnConfig {
            target,
            context_idx,
            subject,
            stage,
            enabled,
        }
    })
}

fn parse_target(raw: &str) -> Option<TraceBlastnTarget> {
    let parts: Vec<&str> = raw
        .split(|c: char| c == ',' || c == ':' || c == ';' || c.is_whitespace())
        .filter(|part| !part.is_empty())
        .collect();
    if parts.len() != 4 {
        eprintln!(
            "[TRACE_BLASTN] invalid LOSAT_TRACE_BLASTN_HSP (expected 4 integers): {:?}",
            raw
        );
        return None;
    }
    Some(TraceBlastnTarget {
        q_start: parts[0].parse().ok()?,
        q_end: parts[1].parse().ok()?,
        s_start: parts[2].parse().ok()?,
        s_end: parts[3].parse().ok()?,
    })
}

fn parse_stage(raw: &str) -> Option<StageFilter> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "all" => Some(StageFilter::All),
        "seed" | "lookup" => Some(StageFilter::Seed),
        "ungapped" => Some(StageFilter::Ungapped),
        "prelim" | "preliminary" => Some(StageFilter::Prelim),
        "traceback" => Some(StageFilter::Traceback),
        "purge" | "reeval" | "tree" => Some(StageFilter::Purge),
        "hitlist" | "final" => Some(StageFilter::Hitlist),
        other => {
            eprintln!(
                "[TRACE_BLASTN] invalid LOSAT_TRACE_BLASTN_STAGE {:?}; using all",
                other
            );
            None
        }
    }
}

#[inline]
pub fn enabled() -> bool {
    config().enabled
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3918
// ```c
// tmp_hsp.context = context;
// tmp_hsp.query.offset = q_start;
// tmp_hsp.query.end = q_end;
// tmp_hsp.subject.offset = s_start;
// tmp_hsp.subject.end = s_end;
// if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
//                                   hit_options->min_diag_separation))
// ```
pub fn should_trace_range(
    stage: &str,
    context_idx: u32,
    subject_idx: usize,
    subject_id: &str,
    q_offset: usize,
    q_end: usize,
    s_offset: usize,
    s_end: usize,
    query_length: usize,
    query_frame: i32,
) -> bool {
    let cfg = config();
    if !cfg.enabled || !stage_matches(cfg.stage, stage) {
        return false;
    }
    if let Some(wanted) = cfg.context_idx {
        if wanted != context_idx {
            return false;
        }
    }
    if let Some(subject) = cfg.subject.as_deref() {
        if subject != subject_id && subject.parse::<usize>().ok() != Some(subject_idx) {
            return false;
        }
    }
    match cfg.target {
        None => true,
        Some(target) => adjusted_range_overlaps_target(
            target,
            q_offset,
            q_end,
            s_offset,
            s_end,
            query_length,
            query_frame,
        ),
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692
// ```c
// Int4 extra_start =
//     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
// ...
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// Blast_HSPListSortByScore(hsp_list);
// ```
pub fn should_trace_subject(
    stage: &str,
    context_idx: Option<u32>,
    subject_idx: usize,
    subject_id: &str,
) -> bool {
    let cfg = config();
    if !cfg.enabled || !stage_matches(cfg.stage, stage) {
        return false;
    }
    if let (Some(wanted), Some(actual)) = (cfg.context_idx, context_idx) {
        if wanted != actual {
            return false;
        }
    } else if cfg.context_idx.is_some() && context_idx.is_none() {
        return false;
    }
    if let Some(subject) = cfg.subject.as_deref() {
        if subject != subject_id && subject.parse::<usize>().ok() != Some(subject_idx) {
            return false;
        }
    }
    true
}

fn stage_matches(filter: StageFilter, stage: &str) -> bool {
    if filter == StageFilter::All {
        return true;
    }
    matches!(
        (filter, stage),
        (StageFilter::Seed, "seed" | "lookup")
            | (StageFilter::Ungapped, "ungapped")
            | (StageFilter::Prelim, "prelim" | "preliminary")
            | (StageFilter::Traceback, "traceback")
            | (StageFilter::Purge, "purge" | "reeval" | "tree")
            | (StageFilter::Hitlist, "hitlist" | "final")
    )
}

fn adjusted_range_overlaps_target(
    target: TraceBlastnTarget,
    q_offset: usize,
    q_end: usize,
    s_offset: usize,
    s_end: usize,
    query_length: usize,
    query_frame: i32,
) -> bool {
    let (q_start_out, q_end_out, s_start_out, s_end_out) =
        adjusted_offsets(q_offset, q_end, s_offset, s_end, query_length, query_frame);
    ranges_overlap(q_start_out, q_end_out, target.q_start, target.q_end)
        && ranges_overlap(s_start_out, s_end_out, target.s_start, target.s_end)
}

fn adjusted_offsets(
    query_offset: usize,
    query_end: usize,
    subject_offset: usize,
    subject_end: usize,
    query_length: usize,
    query_frame: i32,
) -> (i64, i64, i64, i64) {
    if query_frame < 0 {
        let q_end = query_length.saturating_sub(query_offset) as i64;
        let q_start = query_length.saturating_sub(query_end).saturating_add(1) as i64;
        let s_end = subject_offset.saturating_add(1) as i64;
        let s_start = subject_end as i64;
        (q_start, q_end, s_start, s_end)
    } else {
        (
            query_offset.saturating_add(1) as i64,
            query_end as i64,
            subject_offset.saturating_add(1) as i64,
            subject_end as i64,
        )
    }
}

#[inline]
fn ranges_overlap(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> bool {
    let (a_lo, a_hi) = if a_start <= a_end {
        (a_start, a_end)
    } else {
        (a_end, a_start)
    };
    let (b_lo, b_hi) = if b_start <= b_end {
        (b_start, b_end)
    } else {
        (b_end, b_start)
    };
    a_lo <= b_hi && b_lo <= a_hi
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:436-472
// ```c
// BlastGetOffsetsForGappedAlignment(..., &q_start, &s_start);
// ...
// BlastGetStartForGappedAlignmentNucl(query, subject, hsp);
// ...
// AdjustSubjectRange(&s_start, &adjusted_s_length, q_start,
//                    query_length, &start_shift);
// ```
pub fn log(stage: &str, message: impl AsRef<str>) {
    if enabled() {
        eprintln!("[TRACE_BLASTN] stage={} {}", stage, message.as_ref());
    }
}
