//! EXPERIMENT (LOSAT_X_TBXBUCKET / LOSAT_X_TBXBUCKETSHADOW): the two-hit stage of a unit of the
//! unit path (x_pair_par.rs) with its hits grouped by diagonal range. A child module of
//! x_pair_par.rs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-516
//! ```c
//! while (scan_range[1] <= scan_range[2]) {
//!     /* scan the subject sequence for hits */
//!     hits = scansub(lookup_wrap, subject,
//!                               offset_pairs, array_size, scan_range);
//! ...
//!     for (i = 0; i < hits; ++i) {
//!         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//!         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//! ...
//!         diag_coord = (query_offset - subject_offset) & diag_mask;
//! ```
//! A unit (`scan_unit`) applies the per-hit body (`x_hit`) to its hits in scan order on a fresh
//! unit table. The table of a unit spans the query group and the subject chunk (4 bytes per cell:
//! 1.65 MiB for AP027131 x AP027133, 13 MiB for Sakai x MG1655), and the hits of one scan call are
//! spread over all of it, so every hit is a random access into a table that is larger than the
//! first-level cache, and from about 2 MiB larger than the second.
//!
//! This path keeps the order of the hits per cell and changes the order across cells, as
//! LOSAT_X_SEEDBUCKET (x_seed_bucket.rs) does for the reference scan loop: the hits are appended,
//! in scan order, to one segment per range of `2^shift` consecutive unit cells (16 KiB of table
//! for the default 4,096 cells), and when a segment is full, or the unit's scan has ended, the
//! segments are processed range by range, each in scan order (a flush). The exactness argument is
//! the one at the head of x_seed_bucket.rs:
//!
//! 1. `x_hit` reads and writes only the cell of its hit; the extension it starts reads only the two
//!    sequences and the matrix; the read-only values of the unit (`XHitCtx`) are fixed. So the
//!    states a cell goes through, the extensions it starts and the init HSPs they save depend only
//!    on the hits of that cell, taken in scan order.
//! 2. A cell lies in one range, and a segment holds its range's hits in the order they were
//!    appended, which is scan order; a flush processes every hit appended since the previous flush
//!    before any later hit is appended. So every cell sees its hits in scan order.
//! 3. An init HSP is saved when the hit that starts its extension is processed, so a flush saves
//!    its HSPs range by range. Each appended hit carries its position among the hits appended since
//!    the previous flush (`seq`), the flush records the position of every hit that saves an HSP,
//!    and `restore_flush_order` (x_seed_bucket.rs) puts the HSPs that flush saved back into the
//!    order of the positions. All hits of a flush come after the hits of the flush before it, so
//!    the unit's list is the list of `scan_unit` element for element, before the stable score sort
//!    of the merge.
//!
//! Hit records are 8 bytes: the subject offset (32 bits), then `seq` and the cell's offset within
//! its range (`shift` bits). The query offset is the cell plus the unit's lowest diagonal plus the
//! subject offset (the inverse of the cell computation of `scan_unit`, in wrapping `u32`). A
//! segment has `2^cap_shift` records and there are `ranges` segments, chosen so that
//! `ranges * 2^cap_shift <= 2^(32 - shift)`: a flush holds at most that many hits, so `seq` fits
//! its field.
//!
//! `LOSAT_X_TBXBUCKETSHADOW=1` runs `scan_unit` and this path on every unit (each on a fresh table)
//! and stops at the first difference of the saved records or of the final table.

use super::*;

// No NCBI counterpart: counters of LOSAT_X_TBXBUCKET (units bucketed, flushes, hits) and of its
// shadow; they do not change any value NCBI computes.
static TB_UNITS: AtomicU64 = AtomicU64::new(0);
static TB_FLUSHES: AtomicU64 = AtomicU64::new(0);
static TB_HITS: AtomicU64 = AtomicU64::new(0);
static TB_SHADOW_UNITS: AtomicU64 = AtomicU64::new(0);
static TB_SHADOW_RECORDS: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: which order the units of the unit path use (scheduling only).
/// 0 = scan order (`scan_unit`), 1 = LOSAT_X_TBXBUCKET, 2 = LOSAT_X_TBXBUCKETSHADOW (both,
/// compared). The path needs LOSAT_X_SEEDCOMPACT (the records of `run_unit_compact`).
pub(super) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_TBXBUCKETSHADOW").is_some() {
            2
        } else if std::env::var("LOSAT_X_TBXBUCKET").is_ok_and(|v| v != "0") {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: tuning values of the grouping (prototype; they change speed only).
fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

/// log2 of the cells of a range: 4,096 cells, 16 KiB of table. A constant, so that the record's
/// fields are cut with immediate shifts (sweep 1, 2026-10-18: 2^11 / 2^12 / 2^13 alike).
const SHIFT: u32 = 12;
/// The cell's offset within its range, in a record's low word.
const MASK: u32 = (1 << SHIFT) - 1;

/// Smallest unit table, in cells, that is bucketed (LOSAT_X_TBXBUCKET_MIN_CELLS, default 0).
fn min_cells() -> usize {
    static V: OnceLock<usize> = OnceLock::new();
    *V.get_or_init(|| env_usize("LOSAT_X_TBXBUCKET_MIN_CELLS", 0))
}

/// Buffered hits per 16 cells of the unit table (LOSAT_X_TBXBUCKET_DENSITY, default 4: a flush
/// every quarter of a hit per cell).
fn density16() -> usize {
    static V: OnceLock<usize> = OnceLock::new();
    *V.get_or_init(|| env_usize("LOSAT_X_TBXBUCKET_DENSITY", 4).max(1))
}

/// Fewest buffered hits of a unit (LOSAT_X_TBXBUCKET_MIN_BUDGET, default 16,384).
fn min_budget() -> usize {
    static V: OnceLock<usize> = OnceLock::new();
    *V.get_or_init(|| env_usize("LOSAT_X_TBXBUCKET_MIN_BUDGET", 1 << 14).max(1))
}

/// Per-thread buffers of the grouping, kept across units.
#[derive(Default)]
pub(super) struct XBucketScratch {
    /// `ranges` segments of `2^cap_shift` records each.
    recs: Vec<u64>,
    /// Per segment, the index in `recs` of its next record (segment `b` starts at
    /// `b << cap_shift`).
    w: Vec<u32>,
    /// `seq` of every hit of the current flush that saved an init HSP, in saving order.
    keys: Vec<u32>,
}

/// The segment layout of a unit of `width` cells: (ranges, cap_shift), or `None` when the records
/// cannot hold a flush of at least one hit per segment.
fn layout(width: usize, shift: u32) -> Option<(usize, u32)> {
    let ranges = width.div_ceil(1usize << shift).max(1);
    let budget = (width.saturating_mul(density16()) / 16)
        .max(min_budget())
        .min(1 << 32);
    let per = budget.div_ceil(ranges).max(16).next_power_of_two();
    let seq_limit = 1usize << (32 - shift);
    let mut cap_shift = per.trailing_zeros();
    while cap_shift > 0 && ranges << cap_shift > seq_limit {
        cap_shift -= 1;
    }
    (ranges << cap_shift <= seq_limit).then_some((ranges, cap_shift))
}

/// The summary line (with LOSAT_X_STATS or the shadow).
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub(super) fn print_summary() {
    if mode() == 2 {
        eprintln!(
            "[X_SHADOW] TBXBUCKET units_compared={} records_compared={} (records and unit tables identical)",
            TB_SHADOW_UNITS.load(AtomicOrdering::Relaxed),
            TB_SHADOW_RECORDS.load(AtomicOrdering::Relaxed)
        );
    }
    if mode() != 0 && (std::env::var_os("LOSAT_X_STATS").is_some() || timing()) {
        eprintln!(
            "[X_STATS] TBXBUCKET units={} flushes={} hits={} shift={} density16={} min_cells={}",
            TB_UNITS.load(AtomicOrdering::Relaxed),
            TB_FLUSHES.load(AtomicOrdering::Relaxed),
            TB_HITS.load(AtomicOrdering::Relaxed),
            SHIFT,
            density16(),
            min_cells()
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
// ```c
// for (i = 0; i < hits; ++i) {
//     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
//     diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// The per-hit body applied to the buffered hits range by range, each range in scan order (module
// doc, items 1-3), and the saved init HSPs put back into scan order.
/// Processes and empties every segment, then restores the scan order of the init HSPs saved.
/// Its own function, so that the loop that appends the hits keeps its values in registers.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn flush<S: XInitSink, P: XProbe>(
    recs: &[u64],
    w: &mut [u32],
    keys: &mut Vec<u32>,
    cap_shift: u32,
    lo: u32,
    diag: &mut [DiagStruct],
    h: &XHitCtx<'_>,
    init_hsps: &mut S,
    probe: &mut P,
) {
    TB_FLUSHES.fetch_add(1, AtomicOrdering::Relaxed);
    let start = init_hsps.saved_len();
    let window = h.window;
    let diag_ptr = diag.as_mut_ptr();
    for (b, w_b) in w.iter_mut().enumerate() {
        let begin = b << cap_shift;
        let end = *w_b as usize;
        if end == begin {
            continue;
        }
        let base = (b << SHIFT) as u32;
        for &rec in &recs[begin..end] {
            let subject_offset = (rec >> 32) as u32;
            let low = rec as u32;
            let cell = base | (low & MASK);
            let query_offset = cell.wrapping_add(lo).wrapping_add(subject_offset);
            debug_assert!((cell as usize) < diag.len());
            // SAFETY: the cell was computed by `scan_unit_bucketed` from a hit of the unit, so it
            // is below the unit's width = diag.len().
            let diag_entry = unsafe { &mut *diag_ptr.add(cell as usize) };
            probe.at(cell as usize, subject_offset);
            let outcome = x_hit(
                diag_entry,
                query_offset,
                subject_offset,
                h,
                init_hsps,
                probe,
            );
            if outcome == XHit::Saved {
                keys.push(low >> SHIFT);
            }
            probe.hit(outcome, window);
        }
        *w_b = begin as u32;
    }
    init_hsps.restore_flush_order(start, keys);
    probe.batch_end();
}

impl XPairPar {
    // EXPERIMENT (LOSAT_X_TBXBUCKET): `run_unit_compact` with the bucketed scan of a unit (module
    // doc) for a unit table of at least `min_cells()` cells; smaller tables take `scan_unit`.
    // LOSAT_X_TBXBUCKETSHADOW runs both on every unit and compares. LOSAT_X_SEEDUNMAP acts as in
    // `run_unit_compact`.
    /// The saved init HSPs of one unit as compact records, in its scan order.
    #[inline(never)]
    pub(super) fn run_unit_compact_bucketed(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<XUnitHsp> {
        self.run_unit_compact_bucketed_with(unit, inp, scratch, &mut XNoProbe)
    }

    // Diagnostic (LOSAT_TIMING, LOSAT_X_STATS): `run_unit_compact_bucketed` with a probe.
    /// `run_unit_compact_bucketed`, recording into the thread's probe `P`.
    #[inline(never)]
    pub(super) fn run_unit_compact_bucketed_probed<P: XProbeSlot>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
    ) -> Vec<XUnitHsp> {
        let mut probe = std::mem::take(P::slot(scratch));
        let exact = self.run_unit_compact_bucketed_with(unit, inp, scratch, &mut probe);
        *P::slot(scratch) = probe;
        exact
    }

    #[inline(always)]
    fn run_unit_compact_bucketed_with<P: XProbe>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        probe: &mut P,
    ) -> Vec<XUnitHsp> {
        let mut records = std::mem::take(&mut scratch.records);
        records.clear();
        let mapped = self.x_unit_records(unit, scratch, &mut records);
        if mode() == 2 {
            let mut reference: Vec<XUnitHsp> = Vec::new();
            self.scan_unit(unit, inp, scratch, &mut reference, &mut XNoProbe);
            let reference_diag = scratch.diag.clone();
            self.scan_unit_bucketed(unit, inp, scratch, &mut records, probe);
            if records != reference || scratch.diag != reference_diag {
                let first = records
                    .iter()
                    .zip(&reference)
                    .position(|(a, b)| a != b)
                    .unwrap_or(records.len().min(reference.len()));
                let cell = scratch
                    .diag
                    .iter()
                    .zip(&reference_diag)
                    .position(|(a, b)| a != b);
                panic!(
                    "LOSAT_X_TBXBUCKETSHADOW: subject {} frame {} chunk {} group {}: records {} vs {} (first difference at {}), first table difference at {:?}",
                    inp.s_idx, unit.f, unit.k, unit.g, records.len(), reference.len(), first, cell
                );
            }
            TB_SHADOW_UNITS.fetch_add(1, AtomicOrdering::Relaxed);
            TB_SHADOW_RECORDS.fetch_add(records.len() as u64, AtomicOrdering::Relaxed);
        } else if x_unit_width(&self.groups[unit.g], &unit.chunk) >= min_cells() {
            self.scan_unit_bucketed(unit, inp, scratch, &mut records, probe);
        } else {
            self.scan_unit(unit, inp, scratch, &mut records, probe);
        }
        let exact = x_unit_copy(&records, mapped);
        scratch.records = records;
        exact
    }

    // The scan loop of `scan_unit` (RI:3082-3168 on one unit, its own module doc item 3) with the
    // hits appended to the segments of their ranges instead of processed at once (module doc).
    /// Scans one unit and saves its init HSPs, in scan order, into `init_hsps`.
    #[inline(always)]
    fn scan_unit_bucketed<S: XInitSink, P: XProbe>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
        probe: &mut P,
    ) {
        let group = &self.groups[unit.g];
        let sub = &group.lookup;
        let wordsize = inp.wordsize;
        let window = inp.window;
        let s_frame = &inp.s_frames_preliminary[unit.f];
        let subject_full = &s_frame.aa_seq;
        let subject_all = &subject_full[1..subject_full.len() - 1];
        let chunk = unit.chunk;
        let chunk_end = chunk.offset.saturating_add(chunk.length);
        let subject = &subject_all[chunk.offset..chunk_end];
        if subject.len() < wordsize as usize {
            return;
        }
        let width = x_unit_width(group, &chunk);
        let Some((ranges, cap_shift)) = layout(width, SHIFT) else {
            self.scan_unit(unit, inp, scratch, init_hsps, probe);
            return;
        };

        let offset_array_size: i32 = OFFSET_ARRAY_SIZE + sub.longest_chain.max(0);
        let offset_pairs_len =
            crate::algorithm::tblastx::x_scan_bf::x_offset_pairs_len(offset_array_size);
        if scratch.offset_pairs.len() < offset_pairs_len {
            scratch
                .offset_pairs
                .resize(offset_pairs_len, OffsetPair::default());
        }
        let lo = group.q_min.wrapping_sub(chunk.length as u32 - 1);
        probe.unit_begin(width);
        if scratch.diag.capacity() < width
            && x_unit_mapped(scratch.unmap, x_unit_span(group), chunk.length)
        {
            scratch.diag = Vec::with_capacity(x_mapped_len::<DiagStruct>(width));
            SU_DIAG.fetch_add(1, AtomicOrdering::Relaxed);
        }
        scratch.diag.clear();
        scratch.diag.resize(width, DiagStruct::default());

        let XBucketScratch { recs, w, keys } = &mut scratch.bucket;
        // `ranges << cap_shift <= 2^(32 - SHIFT)` (`layout`): the record indices fit a `u32`.
        let rec_len = ranges << cap_shift;
        if recs.len() < rec_len {
            recs.resize(rec_len, 0);
        }
        w.clear();
        w.extend((0..ranges).map(|b| (b << cap_shift) as u32));
        keys.clear();
        let cap_mask = (1u32 << cap_shift) - 1;
        TB_UNITS.fetch_add(1, AtomicOrdering::Relaxed);

        let h = XHitCtx {
            diag_offset: window,
            window,
            wordsize,
            lookup: sub,
            contexts: inp.contexts,
            cutoff_scores: inp.cutoff_scores,
            x_dropoff_per_context: inp.x_dropoff_per_context,
            subject,
            s_f_idx: unit.f,
            s_frame: s_frame.frame,
            s_idx: inp.s_idx,
            s_len: inp.s_len,
        };

        // Hits appended since the last flush (`seq` of the next one).
        let mut total: u32 = 0;
        let mut hits_sum: u64 = 0;
        let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
        let scan_interiors = tblastx_scan_interiors(chunk.length, None);
        for (interior_start, interior_end) in scan_interiors {
            let seq_ranges = clip_tblastx_seq_ranges_for_scan_interior(
                &base_seq_ranges,
                interior_start,
                interior_end,
                wordsize as usize,
                subject.len(),
            );
            if seq_ranges.is_empty() {
                continue;
            }
            let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];
            while scan_range[1] <= scan_range[2] {
                let prev_scan_left = scan_range[1];
                let hits = crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
                    group.scan_bf,
                    sub,
                    subject,
                    &seq_ranges,
                    &mut scratch.offset_pairs[..],
                    offset_array_size,
                    &mut scan_range,
                );
                if hits == 0 && scan_range[1] == prev_scan_left {
                    break;
                }
                probe.batch(hits as usize);
                hits_sum += hits as u64;
                let pairs = &scratch.offset_pairs[..hits as usize];
                let mut recs_ptr = recs.as_mut_ptr();
                let mut w_ptr = w.as_mut_ptr();
                for pair in pairs {
                    let subject_offset = pair.s_off;
                    let cell = pair.q_off.wrapping_sub(subject_offset).wrapping_sub(lo);
                    debug_assert!((cell as usize) < width, "unit diagonal outside its table");
                    let b = (cell >> SHIFT) as usize;
                    let rec = (u64::from(subject_offset) << 32)
                        | u64::from((total << SHIFT) | (cell & MASK));
                    total += 1;
                    // SAFETY: cell < width, so b < ranges == w.len(); w[b] is below the end of
                    // segment b (a full segment is flushed at once, below), so below
                    // ranges << cap_shift <= recs.len(). The pointers are taken again after
                    // every flush, which borrows `recs` and `w` as slices.
                    let next = unsafe {
                        let idx = *w_ptr.add(b);
                        *recs_ptr.add(idx as usize) = rec;
                        *w_ptr.add(b) = idx + 1;
                        idx + 1
                    };
                    if next & cap_mask == 0 {
                        // Segment b is full.
                        flush(
                            recs,
                            w,
                            keys,
                            cap_shift,
                            lo,
                            &mut scratch.diag,
                            &h,
                            init_hsps,
                            probe,
                        );
                        total = 0;
                        recs_ptr = recs.as_mut_ptr();
                        w_ptr = w.as_mut_ptr();
                    }
                }
            }
        }
        flush(
            recs,
            w,
            keys,
            cap_shift,
            lo,
            &mut scratch.diag,
            &h,
            init_hsps,
            probe,
        );
        TB_HITS.fetch_add(hits_sum, AtomicOrdering::Relaxed);
        probe.unit_end();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // No NCBI counterpart: the segment layout keeps a flush's positions within their field.
    #[test]
    fn layout_keeps_seq_in_its_field() {
        for shift in 8..=16u32 {
            for width in [1usize, 255, 4096, 422_758, 3_380_072, 1 << 26] {
                if let Some((ranges, cap_shift)) = layout(width, shift) {
                    assert!(ranges << shift >= width);
                    assert!(ranges << cap_shift <= 1usize << (32 - shift));
                }
            }
        }
    }
}
