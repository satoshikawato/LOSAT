//! HSP tracing infrastructure for debugging
//!
//! Enable tracing by setting environment variable:
//! - `LOSAT_TRACE_HSP="qstart,qend,sstart,send"`
//!   Example: `LOSAT_TRACE_HSP="111880,111860,163496,163476"`
//!
//! This traces a single HSP through the pipeline by matching the final outfmt6/7
//! coordinates (qstart,qend,sstart,send).

use std::sync::OnceLock;

use super::chaining::UngappedHit;
use crate::common::Hit;

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
    if !crate::utils::xdiag::ON {
        return None;
    }
    *TRACE_HSP_TARGET.get_or_init(|| {
        let raw = crate::utils::xdiag::var("LOSAT_TRACE_HSP").ok()?;
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1095-1101
// (`s_BlastSegGetTranslatedOffsets`, `Int4` arithmetic)
// ```c
//    if (segment->frame < 0)	{
//       *start = seq_length - CODON_LENGTH*segment->offset + segment->frame;
//       *end = seq_length - CODON_LENGTH*segment->end + segment->frame + 1;
//    } else if (segment->frame > 0)	{
//       *start = CODON_LENGTH*(segment->offset) + segment->frame - 1;
//       *end = CODON_LENGTH*segment->end + segment->frame - 2;
// ```
// Diagnostics: the 1-based nucleotide coordinates that `convert_coords`
// (extension/mod.rs) gives the reports for an HSP, with its ends read as NCBI's `Int4`
// (`UngappedHit::q_end_int4`) and computed in `i64`: an HSP of the one-hit word finder of length 0
// or less ends at or before its start, below 0 near a frame's start, which the `usize` form does
// not take. On an HSP of positive length the values are `convert_coords`'s.
/// `(q_start, q_end, s_start, s_end)` of `convert_coords` for `h`, its ends as `Int4`.
pub(crate) fn int4_nucleotide_coords(h: &UngappedHit) -> (i64, i64, i64, i64) {
    fn convert(aa_start: u32, aa_end: i32, frame: i8, dna_len: u32) -> (i64, i64) {
        let (start, end, len) = (i64::from(aa_start), i64::from(aa_end), i64::from(dna_len));
        let shift = i64::from(frame.unsigned_abs()) - 1;
        if frame > 0 {
            (start * 3 + shift + 1, end * 3 + shift)
        } else {
            (len - start * 3 - shift, len - end * 3 - shift + 1)
        }
    }
    let (q_start, q_end) = convert(h.q_aa_start, h.q_end_int4(), h.q_frame, h.q_orig_len);
    let (s_start, s_end) = convert(h.s_aa_start, h.s_end_int4(), h.s_frame, h.s_orig_len);
    (q_start, q_end, s_start, s_end)
}

/// Trace an UngappedHit if it matches the trace target (its ends read as `Int4`:
/// `int4_nucleotide_coords`).
#[inline]
pub fn trace_ungapped_hit_if_match(stage: &str, h: &UngappedHit) {
    let Some(target) = trace_hsp_target() else {
        return;
    };

    let (q_start, q_end, s_start, s_end) = int4_nucleotide_coords(h);
    if (q_start, q_end, s_start, s_end)
        == (
            i64::from(target.q_start),
            i64::from(target.q_end),
            i64::from(target.s_start),
            i64::from(target.s_end),
        )
    {
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
            stage, h.raw_score, h.bit_score, h.e_value, h.q_start, h.q_end, h.s_start, h.s_end
        );
    }
}
