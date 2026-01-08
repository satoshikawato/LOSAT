//! HSP tracing infrastructure for debugging
//!
//! Enable tracing by setting environment variable:
//! - `LOSAT_TRACE_HSP="qstart,qend,sstart,send"`
//!   Example: `LOSAT_TRACE_HSP="111880,111860,163496,163476"`
//!
//! This traces a single HSP through the pipeline by matching the final outfmt6/7
//! coordinates (qstart,qend,sstart,send).

use std::sync::OnceLock;

use crate::common::Hit;
use super::chaining::UngappedHit;
use super::extension::convert_coords;

// ---------------------------------------------------------------------------
// HSP tracing (opt-in via env var)
// ---------------------------------------------------------------------------

/// Trace a single HSP through the pipeline by matching the final outfmt6/7
/// coordinates (qstart,qend,sstart,send).
#[derive(Clone, Copy, Debug)]
pub struct TraceHspTarget {
    pub q_start: i32,
    pub q_end: i32,
    pub s_start: i32,
    pub s_end: i32,
}

static TRACE_HSP_TARGET: OnceLock<Option<TraceHspTarget>> = OnceLock::new();

/// Get the trace target from environment variable, if set.
#[inline(always)]
pub fn trace_hsp_target() -> Option<TraceHspTarget> {
    *TRACE_HSP_TARGET.get_or_init(|| {
        let raw = std::env::var("LOSAT_TRACE_HSP").ok()?;
        let parts: Vec<&str> = raw
            .split(|c: char| c == ',' || c == ':' || c == ';' || c.is_whitespace())
            .filter(|p| !p.is_empty())
            .collect();
        if parts.len() != 4 {
            eprintln!(
                "[TRACE_HSP] invalid LOSAT_TRACE_HSP (expected 4 integers): {:?}",
                raw
            );
            return None;
        }
        let q_start: i32 = parts[0].parse().ok()?;
        let q_end: i32 = parts[1].parse().ok()?;
        let s_start: i32 = parts[2].parse().ok()?;
        let s_end: i32 = parts[3].parse().ok()?;
        Some(TraceHspTarget {
            q_start,
            q_end,
            s_start,
            s_end,
        })
    })
}

/// Check if the given coordinates match the trace target.
#[inline(always)]
pub fn trace_match_target(
    target: TraceHspTarget,
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
) -> bool {
    q_start as i32 == target.q_start
        && q_end as i32 == target.q_end
        && s_start as i32 == target.s_start
        && s_end as i32 == target.s_end
}

/// Trace an UngappedHit if it matches the trace target.
#[inline]
pub fn trace_ungapped_hit_if_match(stage: &str, h: &UngappedHit) {
    let Some(target) = trace_hsp_target() else {
        return;
    };

    let (q_start, q_end) = convert_coords(h.q_aa_start, h.q_aa_end, h.q_frame, h.q_orig_len);
    let (s_start, s_end) = convert_coords(h.s_aa_start, h.s_aa_end, h.s_frame, h.s_orig_len);

    if trace_match_target(target, q_start, q_end, s_start, s_end) {
        eprintln!(
            "[TRACE_HSP] stage={} raw_score={} e={} ctx_idx={} s_f_idx={} q_frame={} s_frame={} linked_set={} start_of_chain={} q={}-{} s={}-{}",
            stage,
            h.raw_score,
            h.e_value,
            h.ctx_idx,
            h.s_f_idx,
            h.q_frame,
            h.s_frame,
            h.linked_set,
            h.start_of_chain,
            q_start,
            q_end,
            s_start,
            s_end
        );
    }
}

/// Trace a final Hit if it matches the trace target.
#[inline]
pub fn trace_final_hit_if_match(stage: &str, h: &Hit) {
    let Some(target) = trace_hsp_target() else {
        return;
    };
    if trace_match_target(target, h.q_start, h.q_end, h.s_start, h.s_end) {
        eprintln!(
            "[TRACE_HSP] stage={} raw_score={} bit={:.1} e={} q={}-{} s={}-{}",
            stage,
            h.raw_score,
            h.bit_score,
            h.e_value,
            h.q_start,
            h.q_end,
            h.s_start,
            h.s_end
        );
    }
}
