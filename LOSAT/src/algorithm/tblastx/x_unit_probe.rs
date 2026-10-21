//! Diagnostic (LOSAT_TIMING, LOSAT_X_STATS): what the scan of a unit (`scan_unit` and `x_hit` of
//! x_pair_par.rs) records besides its init HSPs. A child module of x_pair_par.rs.
//!
//! The scan is generic over a probe (`XProbe`), which it calls at the start and end of the unit,
//! after each `scansub` batch, around each hit and around each two-hit extension:
//!
//! - `XNoProbe` records nothing. Every call on it is an empty inline function, so the scan of a
//!   run without the two variables is the scan without probes.
//! - `XTimeProbe` (LOSAT_TIMING) splits the units' time into the seed stage (the table's reset,
//!   the scan, the two-hit test, the saving of init HSPs) and the two-hit extensions (the calls
//!   of `extend_hit_two_hit`), and counts the hits and extensions:
//!
//!   `[TIMING] x_unit: unit_sum ..s seed ..s extend ..s (units=N hits=N extends=N
//!   seed_ns_per_hit .. extend_ns_per_call .. clock_overhead ..s table_kib_mean .. table_kib_max ..)`
//!
//!   `unit_sum` adds the units' own wall times over the threads that ran them (`units` of
//!   `[TIMING] x_pairpar` is the wall time of the subjects' unit runs), and seed + extend =
//!   unit_sum. Each extension is timed by two reads of a cycle counter around the call (the
//!   time-stamp counter on x86-64, `Instant` elsewhere), converted to seconds with the units' own
//!   wall times. The cost of a read is measured once (`read_cost`); one read per extension is
//!   taken out of `extend` and one out of `seed`, and `clock_overhead` is the time taken out.
//!   The split is an estimate: the counter is not serializing.
//! - `XStatsProbe` (LOSAT_X_STATS; it takes the place of the timer when both are set) counts, per
//!   cell of the unit's diagonal table (one cell is one diagonal of the concatenated query; with
//!   LOSAT_X_TBXUNITMASK's masked tables, the diagonals that share the cell):
//!   the hits per cell over the unit and over each `scansub` batch (`h̄` and its histogram,
//!   1 / 2 / 3 / 4-7 / 8+, by cells and by hits), the isolated hits (no other hit on the cell
//!   less than `window` subject positions before or after it), the runs of consecutive
//!   `Masked` and of consecutive `TooFar` outcomes on a cell (1 / 2 / 3 / 4-7 / 8-15 / 16-63 /
//!   64+), the outcomes of `x_hit`, and the cells the extensions score
//!   (`extension/two_hit.rs::x_two_hit_cells`, which repeats each extension) against the distinct
//!   subject positions they cover per cell (`rescan` = cells / distinct; the distinct positions
//!   are counted past the furthest position a previous extension of the cell reached). Lines
//!   `[X_STATS] TBXSEED ...`.
//!
//! No NCBI counterpart: timers and counters; they compute no value NCBI computes and none the
//! search reads. The scan's statements, its order and its saved init HSPs are the same with every
//! probe.

use super::{XHit, XScratch};
use crate::algorithm::tblastx::extension::x_two_hit_cells;
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

/// Which probe the units run with: 0 none, 1 `XTimeProbe`, 2 `XStatsProbe`.
pub(super) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            2
        } else if std::env::var_os("LOSAT_TIMING").is_some() {
            1
        } else {
            0
        }
    })
}

/// The calls the scan of a unit makes on its probe.
pub(super) trait XProbe {
    /// What `ext_begin` hands to `ext_end`.
    type Mark: Copy;
    /// The unit starts; its table has `width` cells.
    fn unit_begin(&mut self, width: usize);
    /// `scansub` returned a batch of `hits` hits.
    fn batch(&mut self, hits: usize);
    /// The next hit is on cell `cell` at subject offset `s`.
    fn at(&mut self, cell: usize, s: u32);
    /// The hit took `outcome`; `window` is the two-hit window.
    fn hit(&mut self, outcome: XHit, window: i32);
    /// The hits of the batch are done.
    fn batch_end(&mut self);
    /// An extension starts.
    fn ext_begin(&mut self) -> Self::Mark;
    /// The extension `extend_hit_two_hit(q, s, s_left, s_right, q_right, x_drop, false)` ended.
    #[allow(clippy::too_many_arguments)]
    fn ext_end(
        &mut self,
        mark: Self::Mark,
        q: &[u8],
        s: &[u8],
        s_left: usize,
        s_right: usize,
        q_right: usize,
        x_drop: i32,
    );
    /// The unit ends.
    fn unit_end(&mut self);
}

/// The probe of a run without LOSAT_TIMING and LOSAT_X_STATS: nothing.
#[derive(Default)]
pub(super) struct XNoProbe;

impl XProbe for XNoProbe {
    type Mark = ();
    #[inline(always)]
    fn unit_begin(&mut self, _width: usize) {}
    #[inline(always)]
    fn batch(&mut self, _hits: usize) {}
    #[inline(always)]
    fn at(&mut self, _cell: usize, _s: u32) {}
    #[inline(always)]
    fn hit(&mut self, _outcome: XHit, _window: i32) {}
    #[inline(always)]
    fn batch_end(&mut self) {}
    #[inline(always)]
    fn ext_begin(&mut self) {}
    #[inline(always)]
    fn ext_end(
        &mut self,
        _mark: (),
        _q: &[u8],
        _s: &[u8],
        _s_left: usize,
        _s_right: usize,
        _q_right: usize,
        _x_drop: i32,
    ) {
    }
    #[inline(always)]
    fn unit_end(&mut self) {}
}

/// Where a thread's scratch keeps a probe between its units.
pub(super) trait XProbeSlot: XProbe + Default {
    fn slot(scratch: &mut XScratch) -> &mut Self;
}

// ---------------------------------------------------------------------------------------------
// LOSAT_TIMING

/// A cycle counter: the time-stamp counter on x86-64, nanoseconds of `Instant` elsewhere.
#[inline(always)]
fn ticks() -> u64 {
    #[cfg(target_arch = "x86_64")]
    {
        #[allow(unused_unsafe)]
        // SAFETY: RDTSC reads a counter; it has no memory effect and every x86-64 CPU has it.
        unsafe {
            core::arch::x86_64::_rdtsc()
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        static BASE: OnceLock<Instant> = OnceLock::new();
        BASE.get_or_init(Instant::now).elapsed().as_nanos() as u64
    }
}

/// The ticks between two reads of `ticks` one after the other (the median of 1,023 pairs).
fn read_cost() -> u64 {
    static COST: OnceLock<u64> = OnceLock::new();
    *COST.get_or_init(|| {
        let mut pairs: Vec<u64> = (0..1023)
            .map(|_| {
                let a = ticks();
                let b = ticks();
                b.wrapping_sub(a)
            })
            .collect();
        pairs.sort_unstable();
        pairs[pairs.len() / 2]
    })
}

#[derive(Default)]
struct XTimeTotals {
    units: u64,
    unit_ns: u64,
    unit_ticks: u64,
    ext_ticks: u64,
    extends: u64,
    hits: u64,
    table_cells: u64,
    table_cells_max: u64,
}

static TIME_TOTALS: Mutex<XTimeTotals> = Mutex::new(XTimeTotals {
    units: 0,
    unit_ns: 0,
    unit_ticks: 0,
    ext_ticks: 0,
    extends: 0,
    hits: 0,
    table_cells: 0,
    table_cells_max: 0,
});

/// LOSAT_TIMING: the clocks of one unit, added to `TIME_TOTALS` at its end.
#[derive(Default)]
pub(super) struct XTimeProbe {
    start: Option<Instant>,
    start_ticks: u64,
    width: usize,
    ext_ticks: u64,
    extends: u64,
    hits: u64,
}

impl XProbe for XTimeProbe {
    type Mark = u64;
    fn unit_begin(&mut self, width: usize) {
        *self = XTimeProbe {
            start: Some(Instant::now()),
            start_ticks: ticks(),
            width,
            ..XTimeProbe::default()
        };
    }
    #[inline(always)]
    fn batch(&mut self, hits: usize) {
        self.hits += hits as u64;
    }
    #[inline(always)]
    fn at(&mut self, _cell: usize, _s: u32) {}
    #[inline(always)]
    fn hit(&mut self, _outcome: XHit, _window: i32) {}
    #[inline(always)]
    fn batch_end(&mut self) {}
    #[inline(always)]
    fn ext_begin(&mut self) -> u64 {
        ticks()
    }
    #[inline(always)]
    fn ext_end(
        &mut self,
        mark: u64,
        _q: &[u8],
        _s: &[u8],
        _s_left: usize,
        _s_right: usize,
        _q_right: usize,
        _x_drop: i32,
    ) {
        self.ext_ticks += ticks().wrapping_sub(mark);
        self.extends += 1;
    }
    fn unit_end(&mut self) {
        let unit_ticks = ticks().wrapping_sub(self.start_ticks);
        let unit_ns = self.start.map_or(0, |t| t.elapsed().as_nanos() as u64);
        let mut t = TIME_TOTALS.lock().unwrap_or_else(|e| e.into_inner());
        t.units += 1;
        t.unit_ns += unit_ns;
        t.unit_ticks += unit_ticks;
        t.ext_ticks += self.ext_ticks;
        t.extends += self.extends;
        t.hits += self.hits;
        t.table_cells += self.width as u64;
        t.table_cells_max = t.table_cells_max.max(self.width as u64);
    }
}

impl XProbeSlot for XTimeProbe {
    fn slot(scratch: &mut XScratch) -> &mut Self {
        &mut scratch.time_probe
    }
}

/// The `[TIMING] x_unit` line, when the units ran with `XTimeProbe`.
pub(super) fn timing_line() -> Option<String> {
    if mode() != 1 {
        return None;
    }
    let t = TIME_TOTALS.lock().unwrap_or_else(|e| e.into_inner());
    let ns_per_tick = if t.unit_ticks > 0 {
        t.unit_ns as f64 / t.unit_ticks as f64
    } else {
        0.0
    };
    let overhead_ns = (t.extends * read_cost()) as f64 * ns_per_tick;
    let extend_ns = (t.ext_ticks as f64 * ns_per_tick - overhead_ns).max(0.0);
    let seed_ns = (t.unit_ns as f64 - extend_ns - 2.0 * overhead_ns).max(0.0);
    let per = |ns: f64, n: u64| if n > 0 { ns / n as f64 } else { 0.0 };
    Some(format!(
        "[TIMING] x_unit: unit_sum {:.3}s seed {:.3}s extend {:.3}s (units={} hits={} extends={} seed_ns_per_hit {:.2} extend_ns_per_call {:.1} clock_overhead {:.3}s table_kib_mean {:.0} table_kib_max {:.0})",
        t.unit_ns as f64 / 1e9,
        seed_ns / 1e9,
        extend_ns / 1e9,
        t.units,
        t.hits,
        t.extends,
        per(seed_ns, t.hits),
        per(extend_ns, t.extends),
        2.0 * overhead_ns / 1e9,
        per(t.table_cells as f64 * 4.0 / 1024.0, t.units),
        t.table_cells_max as f64 * 4.0 / 1024.0
    ))
}

// ---------------------------------------------------------------------------------------------
// LOSAT_X_STATS

/// Histogram of hits per cell: 1, 2, 3, 4-7, 8+.
const H_BUCKETS: [&str; 5] = ["1", "2", "3", "4-7", "8+"];
/// Histogram of run lengths: 1, 2, 3, 4-7, 8-15, 16-63, 64+.
const R_BUCKETS: [&str; 7] = ["1", "2", "3", "4-7", "8-15", "16-63", "64+"];

fn h_bucket(n: u32) -> usize {
    match n {
        0 | 1 => 0,
        2 => 1,
        3 => 2,
        4..=7 => 3,
        _ => 4,
    }
}

fn r_bucket(n: u32) -> usize {
    match n {
        0 | 1 => 0,
        2 => 1,
        3 => 2,
        4..=7 => 3,
        8..=15 => 4,
        16..=63 => 5,
        _ => 6,
    }
}

/// Run kinds: 0 none, 1 `Masked`, 2 `TooFar`.
fn run_kind(outcome: XHit) -> u8 {
    match outcome {
        XHit::Masked => 1,
        XHit::TooFar => 2,
        _ => 0,
    }
}

const OUTCOMES: [(XHit, &str); 7] = [
    (XHit::Masked, "masked"),
    (XHit::FlagReset, "flag_reset"),
    (XHit::TooFar, "too_far"),
    (XHit::Overlap, "overlap"),
    (XHit::ContextBoundary, "context_boundary"),
    (XHit::Extended, "extended"),
    (XHit::Saved, "saved"),
];

fn outcome_index(outcome: XHit) -> usize {
    match outcome {
        XHit::Masked => 0,
        XHit::FlagReset => 1,
        XHit::TooFar => 2,
        XHit::Overlap => 3,
        XHit::ContextBoundary => 4,
        XHit::Extended => 5,
        XHit::Saved => 6,
    }
}

/// The counters of the statistics (one unit's, or the run's totals).
#[derive(Clone, Copy, Default)]
struct XStatTotals {
    units: u64,
    hits: u64,
    table_cells: u64,
    table_cells_max: u64,
    outcomes: [u64; 7],
    unit_cells: u64,
    unit_cell_hist: [u64; 5],
    unit_hit_hist: [u64; 5],
    batches: u64,
    batch_cells: u64,
    batch_cell_hist: [u64; 5],
    batch_hit_hist: [u64; 5],
    isolated: u64,
    /// Runs per kind (`Masked`, `TooFar`) and length bucket.
    runs: [[u64; 7]; 2],
    run_hits: [u64; 2],
    ext_calls: u64,
    ext_cells: u64,
    ext_left: u64,
    ext_right: u64,
    ext_distinct: u64,
}

impl XStatTotals {
    fn add(&mut self, o: &XStatTotals) {
        let add5 = |a: &mut [u64; 5], b: &[u64; 5]| a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
        self.units += o.units;
        self.hits += o.hits;
        self.table_cells += o.table_cells;
        self.table_cells_max = self.table_cells_max.max(o.table_cells_max);
        self.outcomes
            .iter_mut()
            .zip(&o.outcomes)
            .for_each(|(x, y)| *x += y);
        self.unit_cells += o.unit_cells;
        add5(&mut self.unit_cell_hist, &o.unit_cell_hist);
        add5(&mut self.unit_hit_hist, &o.unit_hit_hist);
        self.batches += o.batches;
        self.batch_cells += o.batch_cells;
        add5(&mut self.batch_cell_hist, &o.batch_cell_hist);
        add5(&mut self.batch_hit_hist, &o.batch_hit_hist);
        self.isolated += o.isolated;
        for k in 0..2 {
            self.runs[k]
                .iter_mut()
                .zip(&o.runs[k])
                .for_each(|(x, y)| *x += y);
            self.run_hits[k] += o.run_hits[k];
        }
        self.ext_calls += o.ext_calls;
        self.ext_cells += o.ext_cells;
        self.ext_left += o.ext_left;
        self.ext_right += o.ext_right;
        self.ext_distinct += o.ext_distinct;
    }
}

static STAT_TOTALS: Mutex<Option<XStatTotals>> = Mutex::new(None);

/// What the statistics keep per cell of a unit's table.
#[derive(Clone, Copy, Default)]
struct XCellStat {
    /// Hits on the cell in the unit.
    unit_hits: u32,
    /// Hits on the cell in the current batch.
    batch_hits: u32,
    /// Subject offset + 1 of the cell's last hit (0: none yet).
    last_s1: u32,
    /// Length of the cell's current run.
    run_len: u32,
    /// One past the furthest subject position an extension of the cell scored.
    ext_hi: u32,
    /// Kind of the cell's current run (`run_kind`).
    run_kind: u8,
    /// Whether the cell's last hit has a hit less than `window` before it.
    last_near: bool,
}

/// LOSAT_X_STATS: the per-cell state of one unit and its counters.
#[derive(Default)]
pub(super) struct XStatsProbe {
    cells: Vec<XCellStat>,
    unit_touched: Vec<u32>,
    batch_touched: Vec<u32>,
    cur: usize,
    cur_s: u32,
    t: XStatTotals,
}

impl XStatsProbe {
    fn close_run(t: &mut XStatTotals, c: &mut XCellStat) {
        if c.run_kind != 0 {
            let k = c.run_kind as usize - 1;
            t.runs[k][r_bucket(c.run_len)] += 1;
            t.run_hits[k] += c.run_len as u64;
        }
        c.run_kind = 0;
        c.run_len = 0;
    }
}

impl XProbe for XStatsProbe {
    type Mark = ();
    fn unit_begin(&mut self, width: usize) {
        self.cells.clear();
        self.cells.resize(width, XCellStat::default());
        self.unit_touched.clear();
        self.batch_touched.clear();
        self.t = XStatTotals {
            units: 1,
            table_cells: width as u64,
            table_cells_max: width as u64,
            ..XStatTotals::default()
        };
    }
    fn batch(&mut self, hits: usize) {
        self.t.hits += hits as u64;
        self.t.batches += 1;
    }
    fn at(&mut self, cell: usize, s: u32) {
        self.cur = cell;
        self.cur_s = s;
    }
    fn hit(&mut self, outcome: XHit, window: i32) {
        let (cell, s) = (self.cur, self.cur_s);
        let t = &mut self.t;
        t.outcomes[outcome_index(outcome)] += 1;
        let c = &mut self.cells[cell];
        if c.unit_hits == 0 {
            self.unit_touched.push(cell as u32);
        }
        c.unit_hits += 1;
        if c.batch_hits == 0 {
            self.batch_touched.push(cell as u32);
        }
        c.batch_hits += 1;
        // Isolation: the previous hit of the cell is isolated when it had no near hit before it
        // and this one is not near it either.
        if c.last_s1 != 0 {
            let near = s.wrapping_sub(c.last_s1 - 1) < window as u32;
            if !near && !c.last_near {
                t.isolated += 1;
            }
            c.last_near = near;
        }
        c.last_s1 = s + 1;
        let kind = run_kind(outcome);
        if kind != 0 && kind == c.run_kind {
            c.run_len += 1;
        } else {
            Self::close_run(t, c);
            if kind != 0 {
                c.run_kind = kind;
                c.run_len = 1;
            }
        }
    }
    fn batch_end(&mut self) {
        let t = &mut self.t;
        for &cell in &self.batch_touched {
            let c = &mut self.cells[cell as usize];
            let b = h_bucket(c.batch_hits);
            t.batch_cells += 1;
            t.batch_cell_hist[b] += 1;
            t.batch_hit_hist[b] += c.batch_hits as u64;
            c.batch_hits = 0;
        }
        self.batch_touched.clear();
    }
    fn ext_begin(&mut self) {}
    fn ext_end(
        &mut self,
        _mark: (),
        q: &[u8],
        s: &[u8],
        s_left: usize,
        s_right: usize,
        q_right: usize,
        x_drop: i32,
    ) {
        let cells = x_two_hit_cells(q, s, s_left, s_right, q_right, x_drop);
        let t = &mut self.t;
        t.ext_calls += 1;
        t.ext_left += cells.left as u64;
        t.ext_right += cells.right as u64;
        t.ext_cells += 3 + (cells.left + cells.right) as u64;
        let c = &mut self.cells[self.cur];
        let (lo, hi) = (cells.s_lo as u32, cells.s_hi as u32);
        t.ext_distinct += hi.saturating_sub(lo.max(c.ext_hi)) as u64;
        c.ext_hi = c.ext_hi.max(hi);
    }
    fn unit_end(&mut self) {
        let t = &mut self.t;
        for &cell in &self.unit_touched {
            let c = &mut self.cells[cell as usize];
            let b = h_bucket(c.unit_hits);
            t.unit_cells += 1;
            t.unit_cell_hist[b] += 1;
            t.unit_hit_hist[b] += c.unit_hits as u64;
            if !c.last_near {
                t.isolated += 1;
            }
            Self::close_run(t, c);
        }
        let mut totals = STAT_TOTALS.lock().unwrap_or_else(|e| e.into_inner());
        totals.get_or_insert_with(XStatTotals::default).add(t);
    }
}

impl XProbeSlot for XStatsProbe {
    fn slot(scratch: &mut XScratch) -> &mut Self {
        &mut scratch.stats_probe
    }
}

/// The `[X_STATS] TBXSEED` lines, when the units ran with `XStatsProbe`.
pub(super) fn stats_lines() -> Vec<String> {
    if mode() != 2 {
        return Vec::new();
    }
    let t = STAT_TOTALS
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .unwrap_or_default();
    let ratio = |a: u64, b: u64| if b > 0 { a as f64 / b as f64 } else { 0.0 };
    let hist = |names: &[&str], values: &[u64]| {
        names
            .iter()
            .zip(values)
            .map(|(n, v)| format!("{n}:{v}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    let mut lines = Vec::new();
    let mut outcomes = format!(
        "[X_STATS] TBXSEED units={} hits={} table_cells_mean={:.0} table_cells_max={}",
        t.units,
        t.hits,
        ratio(t.table_cells, t.units),
        t.table_cells_max
    );
    for (o, name) in OUTCOMES {
        outcomes.push_str(&format!(" {name}={}", t.outcomes[outcome_index(o)]));
    }
    lines.push(outcomes);
    lines.push(format!(
        "[X_STATS] TBXSEED unit cells={} hbar={:.3} cells_by_hits={} hits_by_hits={}",
        t.unit_cells,
        ratio(t.hits, t.unit_cells),
        hist(&H_BUCKETS, &t.unit_cell_hist),
        hist(&H_BUCKETS, &t.unit_hit_hist)
    ));
    lines.push(format!(
        "[X_STATS] TBXSEED batch batches={} cells={} hbar={:.3} cells_by_hits={} hits_by_hits={}",
        t.batches,
        t.batch_cells,
        ratio(t.hits, t.batch_cells),
        hist(&H_BUCKETS, &t.batch_cell_hist),
        hist(&H_BUCKETS, &t.batch_hit_hist)
    ));
    lines.push(format!(
        "[X_STATS] TBXSEED isolated={} isolated_frac={:.4} masked_runs={} masked_run_hits={} masked_run_len={} too_far_runs={} too_far_run_hits={} too_far_run_len={}",
        t.isolated,
        ratio(t.isolated, t.hits),
        t.runs[0].iter().sum::<u64>(),
        t.run_hits[0],
        hist(&R_BUCKETS, &t.runs[0]),
        t.runs[1].iter().sum::<u64>(),
        t.run_hits[1],
        hist(&R_BUCKETS, &t.runs[1])
    ));
    lines.push(format!(
        "[X_STATS] TBXSEED extend calls={} cells={} cells_per_call={:.2} left={} right={} distinct={} rescan={:.3}",
        t.ext_calls,
        t.ext_cells,
        ratio(t.ext_cells, t.ext_calls),
        t.ext_left,
        t.ext_right,
        t.ext_distinct,
        ratio(t.ext_cells, t.ext_distinct)
    ));
    lines
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(p: &mut XStatsProbe, width: usize, batches: &[&[(usize, u32, XHit)]]) -> XStatTotals {
        p.unit_begin(width);
        for batch in batches {
            p.batch(batch.len());
            for &(cell, s, outcome) in batch.iter() {
                p.at(cell, s);
                p.hit(outcome, 40);
            }
            p.batch_end();
        }
        // `unit_end` adds the counts of the cells' last hits and open runs to `p.t`.
        p.unit_end();
        p.t
    }

    // No NCBI counterpart (statistics). Hits per cell over the unit and per batch, isolated
    // hits, and runs, on a hand-made sequence of hits.
    #[test]
    fn the_statistics_count_cells_isolation_and_runs() {
        use XHit::*;
        let mut p = XStatsProbe::default();
        // Cell 0: hits at 0, 100, 110, 120. Cell 1: one hit. Cell 2: masked, masked, too far,
        // too far, too far, extended.
        let t = run(
            &mut p,
            4,
            &[
                &[
                    (0, 0, TooFar),
                    (2, 5, Masked),
                    (2, 6, Masked),
                    (1, 7, TooFar),
                ],
                &[
                    (0, 100, TooFar),
                    (2, 200, TooFar),
                    (2, 300, TooFar),
                    (2, 400, TooFar),
                ],
                &[(0, 110, Overlap), (0, 120, Extended), (2, 410, Extended)],
            ],
        );
        assert_eq!(t.hits, 11);
        assert_eq!(t.unit_cells, 3);
        // Unit: cell 0 has 4 hits (4-7), cell 1 one, cell 2 six (4-7).
        assert_eq!(t.unit_cell_hist, [1, 0, 0, 2, 0]);
        assert_eq!(t.unit_hit_hist, [1, 0, 0, 10, 0]);
        // Batches: {0:1, 2:2, 1:1}, {0:1, 2:3}, {0:2, 2:1}.
        assert_eq!(t.batch_cells, 7);
        assert_eq!(t.batch_cell_hist, [4, 2, 1, 0, 0]);
        assert_eq!(t.batch_hit_hist, [4, 4, 3, 0, 0]);
        // Isolated: cell 0 at 0 (100 has 110 near it); cell 1 at 7; cell 2 at 200 and 300 (5
        // and 6 are near each other, and 400 and 410).
        assert_eq!(t.isolated, 4);
        // Runs: cell 0 TooFar x2 then Overlap (one too-far run of 2); cell 1 TooFar x1; cell 2
        // Masked x2, TooFar x3.
        assert_eq!(t.runs[0], [0, 1, 0, 0, 0, 0, 0]);
        assert_eq!(t.runs[1], [1, 1, 1, 0, 0, 0, 0]);
        assert_eq!(t.run_hits, [2, 6]);
    }

    // No NCBI counterpart (statistics).
    #[test]
    fn buckets_cover_their_ranges() {
        assert_eq!(
            [1, 2, 3, 4, 7, 8, 1000].map(h_bucket),
            [0, 1, 2, 3, 3, 4, 4]
        );
        assert_eq!(
            [1, 2, 3, 4, 7, 8, 15, 16, 63, 64].map(r_bucket),
            [0, 1, 2, 3, 3, 4, 4, 5, 5, 6]
        );
    }
}
