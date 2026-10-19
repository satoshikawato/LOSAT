//! EXPERIMENT (R7_PROBE; line R's R7 probe of plan v4, `v4-R/DESIGN_R7.md`; not for merge): page
//! faults of the megablast table build, measured before the fill. No NCBI counterpart and no
//! change of any value: every write stores 0 over the 0 the fresh tables hold, and the extra
//! arrays are dropped before the fill.
//!
//! `R7_PROBE` is a comma-separated list:
//! - `pf`: write 0 to one element of every 4 KiB page of `hashtable` (first touch by a write: one
//!   fault per page, where the fill's read-then-write takes the zero page and then a copy-on-write).
//! - `np`: the same for `next_pos` (the fill writes it without reading it: its faults leave the fill).
//! - `dbl=<MiB>`: a fresh zeroed array of that size, each page read then written (the fill's pattern
//!   on `hashtable`), timed and freed.
//! - `sgl=<MiB>`: the same with each page only written.
//!
//! Lines on stderr: `[R7] <name> pages= ... s=` (and `alloc_s`, `free_s` for the extra arrays).

use std::sync::OnceLock;
use std::time::Instant;

const PAGE: usize = 4096;

#[derive(Default)]
struct Probe {
    pf: bool,
    np: bool,
    dbl_mib: usize,
    sgl_mib: usize,
}

fn probe() -> &'static Probe {
    static P: OnceLock<Probe> = OnceLock::new();
    P.get_or_init(|| {
        let mut p = Probe::default();
        let v = std::env::var("R7_PROBE").unwrap_or_default();
        for tok in v.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            let (k, val) = tok.split_once('=').unwrap_or((tok, ""));
            match k {
                "pf" => p.pf = true,
                "np" => p.np = true,
                "dbl" => p.dbl_mib = val.parse().unwrap_or(0),
                "sgl" => p.sgl_mib = val.parse().unwrap_or(0),
                _ => eprintln!("[R7] unknown R7_PROBE token {tok}"),
            }
        }
        p
    })
}

/// Writes 0 to the first element of each 4 KiB page `a` spans; returns the pages written. The
/// caller guarantees every element is 0 (fresh `vec![0; n]`).
fn touch_pages(a: &mut [u32], read_first: bool) -> (usize, u32) {
    let base = a.as_mut_ptr() as usize;
    let mut i = 0usize;
    let mut pages = 0usize;
    let mut x = 0u32;
    while i < a.len() {
        // SAFETY: i < a.len().
        unsafe {
            if read_first {
                x ^= std::ptr::read_volatile(a.as_ptr().add(i));
            }
            std::ptr::write_volatile(a.as_mut_ptr().add(i), x);
        }
        pages += 1;
        let next = ((base + i * 4) / PAGE + 1) * PAGE;
        i = (next - base) / 4;
    }
    (pages, x)
}

fn extra(name: &str, mib: usize, read_first: bool) {
    let n = mib << 18;
    let t0 = Instant::now();
    let mut v = vec![0u32; n];
    let t1 = Instant::now();
    let (pages, x) = touch_pages(&mut v, read_first);
    let t2 = Instant::now();
    std::hint::black_box(&v);
    drop(v);
    let t3 = Instant::now();
    std::hint::black_box(x);
    eprintln!(
        "[R7] {name} mib={mib} pages={pages} alloc_s={:.6} s={:.6} free_s={:.6}",
        (t1 - t0).as_secs_f64(),
        (t2 - t1).as_secs_f64(),
        (t3 - t2).as_secs_f64()
    );
}

/// Called by `build_mb_lookup` after the tables are allocated and before the fill.
#[inline(never)]
pub(super) fn before_fill(hashtable: &mut [u32], next_pos: &mut [u32]) {
    let p = probe();
    if !(p.pf || p.np || p.dbl_mib > 0 || p.sgl_mib > 0) {
        return;
    }
    if p.dbl_mib > 0 {
        extra("dbl", p.dbl_mib, true);
    }
    if p.sgl_mib > 0 {
        extra("sgl", p.sgl_mib, false);
    }
    if p.pf {
        let t = Instant::now();
        let (pages, _) = touch_pages(hashtable, false);
        eprintln!("[R7] pf pages={pages} s={:.6}", t.elapsed().as_secs_f64());
    }
    if p.np {
        let t = Instant::now();
        let (pages, _) = touch_pages(next_pos, false);
        eprintln!("[R7] np pages={pages} s={:.6}", t.elapsed().as_secs_f64());
    }
}
