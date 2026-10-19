//! EXPERIMENT (R9_PROBE; line R's R9 probe of plan v4, `v4-R/DESIGN_R9.md`; not for merge): counts
//! of the megablast greedy traceback and of the word extension. No value changes: the hooks only
//! read the rows and offsets the search computes anyway. Printed at exit as `[R9] ...` lines.
//!
//! `R9_PROBE` is a comma-separated list:
//! - `gc`: per row `d` of every traceback `blast_greedy_align` (one per direction): the band
//!   width, the diagonals that survive the X-drop test (the C loop's `s_FindFirstMismatch` calls)
//!   and their slide lengths. At the end of a direction the rows `d <= best_dist` are the body and
//!   the rows `d > best_dist` the tail (searched until every diagonal drops). Each traceback call
//!   (both directions) goes to the power-of-two bucket of its alignment length (as R7's `gt`).
//! - `we`: per word extension of the scan, the bases compared on the left and on the right
//!   (matches plus the mismatch that stops the loop).
//! - `sc`: the traceback rows run `x_greedy_row_avx2_sc`, the same kernel with its second pass
//!   kept scalar (one `tzcnt` / `lzcnt` per diagonal instead of the vectorized emulation).

use std::cell::RefCell;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};

fn flags() -> (bool, bool) {
    static F: OnceLock<(bool, bool)> = OnceLock::new();
    *F.get_or_init(|| {
        let v = std::env::var("R9_PROBE").unwrap_or_default();
        let has = |t: &str| v.split(',').any(|x| x.trim() == t);
        (has("gc"), has("we"))
    })
}

/// `sc`: the greedy rows use the kernel with a scalar pass 2 (`x_greedy_row_avx2_sc`).
#[inline]
pub(crate) fn sc_on() -> bool {
    static F: OnceLock<bool> = OnceLock::new();
    *F.get_or_init(|| {
        std::env::var("R9_PROBE")
            .map(|v| v.split(',').any(|x| x.trim() == "sc"))
            .unwrap_or(false)
    })
}

#[inline(never)]
pub(crate) fn gc_on() -> bool {
    flags().0
}

#[inline(never)]
pub(crate) fn we_on() -> bool {
    flags().1
}

/// Slide buckets: 0, 1-3, 4-7, 8-15, 16-31, 32+.
const SB: usize = 6;
fn sbucket(n: i64) -> usize {
    match n {
        i64::MIN..=0 => 0,
        1..=3 => 1,
        4..=7 => 2,
        8..=15 => 3,
        16..=31 => 4,
        _ => 5,
    }
}
/// Power-of-two buckets: 0, 1, 2-3, 4-7, ...
const PB: usize = 24;
fn pbucket(n: u64) -> usize {
    if n == 0 {
        0
    } else {
        ((u64::BITS - n.leading_zeros()) as usize).min(PB - 1)
    }
}
fn plabel(b: usize) -> String {
    if b == 0 {
        "0".into()
    } else {
        format!("{}-{}", 1u64 << (b - 1), (1u64 << b) - 1)
    }
}

#[derive(Clone, Copy, Default)]
struct Part {
    rows: u64,
    band: u64,
    live: u64,
    slide: u64,
    slides: [u64; SB],
}
impl Part {
    fn add(&mut self, o: &Part) {
        self.rows += o.rows;
        self.band += o.band;
        self.live += o.live;
        self.slide += o.slide;
        for i in 0..SB {
            self.slides[i] += o.slides[i];
        }
    }
}

#[derive(Clone, Copy, Default)]
struct Call {
    dirs: u64,
    seed: u64,
    body: Part,
    tail: Part,
    retry: Part,
    retries: u64,
    aborted_rows: u64,
}
impl Call {
    fn add(&mut self, o: &Call) {
        self.dirs += o.dirs;
        self.seed += o.seed;
        self.body.add(&o.body);
        self.tail.add(&o.tail);
        self.retry.add(&o.retry);
        self.retries += o.retries;
        self.aborted_rows += o.aborted_rows;
    }
}

#[derive(Default)]
struct Totals {
    flushed: bool,
    calls: [u64; PB],
    by_len: [Call; PB],
    // by direction: 0 = right (forward), 1 = left (reverse)
    by_dir: [Call; 2],
    body_rows_hist: [u64; PB],
    tail_rows_hist: [u64; PB],
    body_live_hist: [u64; PB],
    tail_live_hist: [u64; PB],
    last_band_hist: [u64; PB],
    unflushed: Call,
}

struct Row {
    band: u32,
    live: u32,
    slide: u64,
    slides: [u32; SB],
}

#[derive(Default)]
struct State {
    starts: Vec<i32>,
    rows: Vec<Row>,
    reverse: bool,
    in_dir: bool,
    seed: u64,
    call: Call,
}

static TOTALS: Mutex<Option<Totals>> = Mutex::new(None);

fn with_totals<R>(f: impl FnOnce(&mut Totals) -> R) -> R {
    let mut g = TOTALS.lock().unwrap();
    f(g.get_or_insert_with(Totals::default))
}

thread_local! {
    static ST: RefCell<State> = RefCell::new(State::default());
}

/// The start of one direction of a traceback greedy alignment; `seed`: the slide of row 0.
pub(crate) fn gc_begin(reverse: bool, seed: i32) {
    ST.with(|s| {
        let s = &mut *s.borrow_mut();
        if s.in_dir {
            // the previous direction returned without `gc_end` (fence)
            s.call.aborted_rows += s.rows.len() as u64;
        }
        s.rows.clear();
        s.reverse = reverse;
        s.in_dir = true;
        s.seed = seed.max(0) as u64;
    });
}

/// Before row `d`: which diagonals of the band survive the X-drop test, from the row `d - 1`
/// (`previous` holds the band plus one cell on each side), as the C loop decides.
pub(crate) fn gc_row_pre(previous: &[i32], band_len: usize, first: i32, xdrop_score: i32) {
    ST.with(|s| {
        let s = &mut *s.borrow_mut();
        s.starts.clear();
        for i in 0..band_len {
            let mut seq2_index = previous[i + 2].max(previous[i + 1]) + 1;
            seq2_index = seq2_index.max(previous[i]);
            let extent = seq2_index + seq2_index + first + i as i32;
            s.starts.push(if seq2_index < 0 || extent < xdrop_score {
                -1
            } else {
                seq2_index
            });
        }
    });
}

/// After row `d`: the slide of every surviving diagonal (its new offset minus its start).
pub(crate) fn gc_row_post(current: &[i32]) {
    ST.with(|s| {
        let s = &mut *s.borrow_mut();
        let mut row = Row {
            band: s.starts.len() as u32,
            live: 0,
            slide: 0,
            slides: [0; SB],
        };
        for (i, &st) in s.starts.iter().enumerate() {
            if st < 0 {
                continue;
            }
            let n = (current[i] - st) as i64;
            row.live += 1;
            row.slide += n.max(0) as u64;
            row.slides[sbucket(n)] += 1;
        }
        s.rows.push(row);
    });
}

/// The end of one direction: `converged` false is a try that `max_dist` stopped (retried with
/// a doubled `max_dist`); otherwise rows `1..=best_dist` are the body and the rest the tail.
pub(crate) fn gc_end(converged: bool, best_dist: i32) {
    ST.with(|s| {
        let s = &mut *s.borrow_mut();
        s.in_dir = false;
        let mut body = Part::default();
        let mut tail = Part::default();
        for (j, r) in s.rows.iter().enumerate() {
            let d = j as i32 + 1;
            let p = if !converged || d <= best_dist {
                &mut body
            } else {
                &mut tail
            };
            p.rows += 1;
            p.band += r.band as u64;
            p.live += r.live as u64;
            p.slide += r.slide;
            for i in 0..SB {
                p.slides[i] += r.slides[i] as u64;
            }
        }
        let dir = s.reverse as usize;
        let mut c = Call::default();
        if !converged {
            c.retry = body;
            c.retries = 1;
        } else {
            c.dirs = 1;
            c.seed = s.seed;
            c.body = body;
            c.tail = tail;
        }
        let rows = &s.rows;
        with_totals(|t| {
            if converged {
                t.body_rows_hist[pbucket(body.rows)] += 1;
                t.tail_rows_hist[pbucket(tail.rows)] += 1;
                for (j, r) in rows.iter().enumerate() {
                    if (j as i32 + 1) <= best_dist {
                        t.body_live_hist[pbucket(r.live as u64)] += 1;
                    } else {
                        t.tail_live_hist[pbucket(r.live as u64)] += 1;
                    }
                }
                if let Some(r) = rows.last() {
                    t.last_band_hist[pbucket(r.band as u64)] += 1;
                }
            }
            t.by_dir[dir].add(&c);
        });
        s.call.add(&c);
    });
}

/// The end of one traceback call that returned an alignment of `len` columns (`None`: none).
pub(crate) fn gc_flush(len: Option<usize>) {
    ST.with(|s| {
        let s = &mut *s.borrow_mut();
        let c = std::mem::take(&mut s.call);
        let b = pbucket(len.unwrap_or(0) as u64);
        with_totals(|t| {
            t.flushed = true;
            t.calls[b] += 1;
            t.by_len[b].add(&c);
        });
    });
}

fn part_line(tag: &str, p: &Part) -> String {
    format!(
        "{tag}_rows={} {tag}_band={} {tag}_live={} {tag}_slide={} {tag}_slides={}",
        p.rows,
        p.band,
        p.live,
        p.slide,
        p.slides
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join("/")
    )
}

fn call_line(c: &Call) -> String {
    format!(
        "dirs={} seed={} {} {} retries={} {} aborted_rows={}",
        c.dirs,
        c.seed,
        part_line("body", &c.body),
        part_line("tail", &c.tail),
        c.retries,
        part_line("retry", &c.retry),
        c.aborted_rows
    )
}

fn hist_line(name: &str, h: &[u64; PB]) -> String {
    let parts: Vec<String> = (0..PB)
        .filter(|&b| h[b] != 0)
        .map(|b| format!("{}:{}", plabel(b), h[b]))
        .collect();
    format!("[R9] gc hist {name} {}", parts.join(" "))
}

// word extension counters
static WE: [AtomicU64; 64] = [const { AtomicU64::new(0) }; 64];
const WE_LEFT: usize = 0; // +bucket of compared bases (6)
const WE_RIGHT: usize = 8; // +bucket (6)
const WE_TOTAL: usize = 16; // +bucket of left+right compared (6)
const WE_LEFT_N: usize = 24;
const WE_LEFT_CMP: usize = 25;
const WE_LEFT_FULL: usize = 26; // stopped by the bound, not by a mismatch
const WE_RIGHT_N: usize = 27;
const WE_RIGHT_CMP: usize = 28;
const WE_RIGHT_FULL: usize = 29;
const WE_PASS: usize = 30;
const WE_LEFT_ALL: usize = 31; // left alone reached ext_to (no right loop)
const WE_LEFT_HIST: usize = 32; // +ext_left matched (0..=31)

/// The left loop of one word extension stopped after `ext_left` matches (`max_ext_left`: its bound).
pub(crate) fn we_left(ext_left: usize, max_ext_left: usize, ext_to: usize) {
    let r = Ordering::Relaxed;
    let cmp = ext_left + (ext_left < max_ext_left) as usize;
    WE[WE_LEFT_N].fetch_add(1, r);
    WE[WE_LEFT_CMP].fetch_add(cmp as u64, r);
    WE[WE_LEFT + sbucket(cmp as i64)].fetch_add(1, r);
    WE[WE_LEFT_HIST + ext_left.min(31)].fetch_add(1, r);
    if ext_left >= max_ext_left {
        WE[WE_LEFT_FULL].fetch_add(1, r);
    }
    if ext_left >= ext_to {
        WE[WE_LEFT_ALL].fetch_add(1, r);
        WE[WE_TOTAL + sbucket(cmp as i64)].fetch_add(1, r);
        WE[WE_PASS].fetch_add(1, r);
    }
}

/// The right loop stopped after `ext_right` matches of its bound `limit`; the left loop had compared
/// `left_cmp` bases.
pub(crate) fn we_right(ext_left: usize, max_ext_left: usize, ext_right: usize, limit: usize) {
    let r = Ordering::Relaxed;
    let cmp = ext_right + (ext_right < limit) as usize;
    let left_cmp = ext_left + (ext_left < max_ext_left) as usize;
    WE[WE_RIGHT_N].fetch_add(1, r);
    WE[WE_RIGHT_CMP].fetch_add(cmp as u64, r);
    WE[WE_RIGHT + sbucket(cmp as i64)].fetch_add(1, r);
    WE[WE_TOTAL + sbucket((cmp + left_cmp) as i64)].fetch_add(1, r);
    if ext_right >= limit {
        WE[WE_RIGHT_FULL].fetch_add(1, r);
        WE[WE_PASS].fetch_add(1, r);
    }
}

/// Prints the counts (called at exit by `x_mb_lookup_print_summary`).
pub(crate) fn print() {
    let (gc, we) = flags();
    if gc {
        let pending = ST.with(|s| std::mem::take(&mut s.borrow_mut().call));
        with_totals(|t| {
            t.unflushed.add(&pending);
            let mut all = Call::default();
            for b in 0..PB {
                if t.calls[b] == 0 {
                    continue;
                }
                all.add(&t.by_len[b]);
                eprintln!(
                    "[R9] gc len={} calls={} {}",
                    plabel(b),
                    t.calls[b],
                    call_line(&t.by_len[b])
                );
            }
            eprintln!(
                "[R9] gc all calls={} {}",
                t.calls.iter().sum::<u64>(),
                call_line(&all)
            );
            eprintln!("[R9] gc dir=right {}", call_line(&t.by_dir[0]));
            eprintln!("[R9] gc dir=left {}", call_line(&t.by_dir[1]));
            eprintln!("[R9] gc unflushed {}", call_line(&t.unflushed));
            eprintln!("{}", hist_line("body_rows_per_dir", &t.body_rows_hist));
            eprintln!("{}", hist_line("tail_rows_per_dir", &t.tail_rows_hist));
            eprintln!("{}", hist_line("body_live_per_row", &t.body_live_hist));
            eprintln!("{}", hist_line("tail_live_per_row", &t.tail_live_hist));
            eprintln!("{}", hist_line("last_band_per_dir", &t.last_band_hist));
        });
    }
    if we {
        let g = |i: usize| WE[i].load(Ordering::Relaxed);
        let b = |o: usize| {
            (0..SB)
                .map(|i| g(o + i).to_string())
                .collect::<Vec<_>>()
                .join("/")
        };
        eprintln!(
            "[R9] we left n={} cmp={} bound={} all={} buckets={}",
            g(WE_LEFT_N),
            g(WE_LEFT_CMP),
            g(WE_LEFT_FULL),
            g(WE_LEFT_ALL),
            b(WE_LEFT)
        );
        eprintln!(
            "[R9] we right n={} cmp={} full={} buckets={}",
            g(WE_RIGHT_N),
            g(WE_RIGHT_CMP),
            g(WE_RIGHT_FULL),
            b(WE_RIGHT)
        );
        eprintln!("[R9] we total pass={} buckets={}", g(WE_PASS), b(WE_TOTAL));
        eprintln!(
            "[R9] we left_matched {}",
            (0..32)
                .map(|i| g(WE_LEFT_HIST + i).to_string())
                .collect::<Vec<_>>()
                .join("/")
        );
    }
}
