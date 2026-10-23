//! PROBE ONLY (v4-B s10, B14 premise check; never merged). A counting global allocator that
//! tags every allocation with the allocating thread and a region tag (set by guards in the
//! code), and keeps live bytes per (thread, tag). Snapshots: at the live-bytes peak (256 KiB
//! steps), at the RssAnon peak (sampler thread, 2 ms, /proc/self/statm), and at phase marks.
//! Report to the file in PROBE_MEM_OUT (or stderr) from `report()`.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicI64, AtomicU64, AtomicU8, AtomicUsize, Ordering::Relaxed};
use std::time::Instant;

pub const NT: usize = 160;
pub const NG: usize = 16;
pub const NPH: usize = 8;

pub const OTHER: u8 = 0;
pub const SCRATCH_NEW: u8 = 1;
pub const PRELIM: u8 = 2;
pub const TRACE: u8 = 3;
pub const TREE: u8 = 4;
pub const DBBLOCK: u8 = 5;
pub const LOOKUP: u8 = 6;
pub const COLLECT: u8 = 7;
pub const REPORT: u8 = 8;
pub const LOOKUP_POOL: u8 = 9;
pub const GREEDY: u8 = 10;
pub const XDROP: u8 = 11;
pub const SCAN: u8 = 12;
pub const PRELIM_OUT: u8 = 13;
pub const SCRATCH_HELD: u8 = 14;
pub const TAG_NAMES: [&str; NG] = [
    "other", "scratch_new", "prelim", "trace", "tree", "dbblock", "lookup", "collect", "report",
    "lookup_pool", "greedy", "xdrop", "diag", "prelim_out", "scratch_held", "query",
];
pub const PHASE_NAMES: [&str; NPH] = [
    "start", "lookup_done", "prelim_start", "prelim_done", "trace_done", "report_start", "end", "p7",
];

#[allow(clippy::declare_interior_mutable_const)]
const ZI: AtomicI64 = AtomicI64::new(0);
#[allow(clippy::declare_interior_mutable_const)]
const ZU: AtomicU64 = AtomicU64::new(0);
#[allow(clippy::declare_interior_mutable_const)]
const ROWI: [AtomicI64; NG] = [ZI; NG];
#[allow(clippy::declare_interior_mutable_const)]
const ROWU: [AtomicU64; NG] = [ZU; NG];
#[allow(clippy::declare_interior_mutable_const)]
const TABI: [[AtomicI64; NG]; NT] = [ROWI; NT];

static LIVE: [[AtomicI64; NG]; NT] = [ROWI; NT];
static CUM: [[AtomicU64; NG]; NT] = [ROWU; NT];
static CNT: [[AtomicU64; NG]; NT] = [ROWU; NT];
static TOTAL: AtomicI64 = AtomicI64::new(0);
static PEAK: AtomicI64 = AtomicI64::new(0);
static SNAP_AT: AtomicI64 = AtomicI64::new(0);
static SNAPPING: AtomicBool = AtomicBool::new(false);
static SNAP: [[AtomicI64; NG]; NT] = [ROWI; NT];
static SNAP_INFO: [AtomicU64; 8] = [ZU; 8]; // total, phase, t_us, mi.arena, mi.uord, mi.ford, mi.hblkhd, anon_kb
static RSS_SNAP: [[AtomicI64; NG]; NT] = [ROWI; NT];
static RSS_INFO: [AtomicU64; 8] = [ZU; 8]; // anon_kb, phase, t_us, live total, resident_kb, mi.arena, mi.uord, mi.ford
static PHASE_SNAP: [[[AtomicI64; NG]; NT]; NPH] = [TABI; NPH];
#[allow(clippy::declare_interior_mutable_const)]
const ROW4: [AtomicU64; 4] = [ZU; 4];
static PHASE_INFO: [[AtomicU64; 4]; NPH] = [ROW4; NPH]; // set, t_us, live total, anon_kb
static PHASE: AtomicU8 = AtomicU8::new(0);
static NEXT_TID: AtomicUsize = AtomicUsize::new(0);
static T0_NS: AtomicU64 = AtomicU64::new(0);
const STEP: i64 = 256 << 10;

thread_local! {
    static TAG: Cell<u8> = const { Cell::new(0) };
    static TID: Cell<u16> = const { Cell::new(u16::MAX) };
}

#[repr(C)]
#[derive(Default, Clone, Copy)]
struct Mallinfo2 {
    arena: usize,
    ordblks: usize,
    smblks: usize,
    hblks: usize,
    hblkhd: usize,
    usmblks: usize,
    fsmblks: usize,
    uordblks: usize,
    fordblks: usize,
    keepcost: usize,
}
extern "C" {
    fn mallinfo2() -> Mallinfo2;
    fn open(path: *const u8, flags: i32, ...) -> i32;
    fn read(fd: i32, buf: *mut u8, n: usize) -> isize;
    fn close(fd: i32) -> i32;
    fn clock_gettime(clk: i32, ts: *mut [i64; 2]) -> i32;
}

fn now_us() -> u64 {
    let mut ts = [0i64; 2];
    unsafe { clock_gettime(1, &mut ts) };
    let ns = ts[0] as u64 * 1_000_000_000 + ts[1] as u64;
    ns.saturating_sub(T0_NS.load(Relaxed)) / 1000
}

/// (resident_kb, anon_kb) from /proc/self/statm without allocating.
fn statm() -> (u64, u64) {
    let mut buf = [0u8; 256];
    let fd = unsafe { open(b"/proc/self/statm\0".as_ptr(), 0) };
    if fd < 0 {
        return (0, 0);
    }
    let n = unsafe { read(fd, buf.as_mut_ptr(), buf.len()) };
    unsafe { close(fd) };
    if n <= 0 {
        return (0, 0);
    }
    let mut vals = [0u64; 3];
    let mut i = 0;
    let mut k = 0;
    while (i as isize) < n && k < 3 {
        let c = buf[i];
        if c.is_ascii_digit() {
            vals[k] = vals[k] * 10 + (c - b'0') as u64;
        } else {
            k += 1;
        }
        i += 1;
    }
    let page_kb = 4;
    (vals[1] * page_kb, (vals[1].saturating_sub(vals[2])) * page_kb)
}

fn tid() -> usize {
    TID.try_with(|t| {
        let v = t.get();
        if v != u16::MAX {
            v as usize
        } else {
            let n = NEXT_TID.fetch_add(1, Relaxed).min(NT - 1);
            t.set(n as u16);
            n
        }
    })
    .unwrap_or(NT - 1)
}

fn cur_tag() -> u8 {
    TAG.try_with(|t| t.get()).unwrap_or(0)
}

pub struct TagGuard(u8);
impl Drop for TagGuard {
    fn drop(&mut self) {
        let _ = TAG.try_with(|t| t.set(self.0));
    }
}
#[inline(never)]
pub fn tag(g: u8) -> TagGuard {
    let old = TAG.try_with(|t| t.replace(g)).unwrap_or(0);
    TagGuard(old)
}

fn copy_tab(dst: &[[AtomicI64; NG]; NT]) {
    for t in 0..NT {
        for g in 0..NG {
            dst[t][g].store(LIVE[t][g].load(Relaxed), Relaxed);
        }
    }
}

#[inline(never)]
fn snapshot(tot: i64) {
    if SNAPPING.swap(true, Relaxed) {
        return;
    }
    SNAP_AT.store(tot + STEP, Relaxed);
    copy_tab(&SNAP);
    let mi = unsafe { mallinfo2() };
    let (_, anon) = statm();
    let v = [
        tot as u64,
        PHASE.load(Relaxed) as u64,
        now_us(),
        mi.arena as u64,
        mi.uordblks as u64,
        mi.fordblks as u64,
        mi.hblkhd as u64,
        anon,
    ];
    for (a, b) in SNAP_INFO.iter().zip(v) {
        a.store(b, Relaxed);
    }
    SNAPPING.store(false, Relaxed);
}

#[inline]
fn add(t: usize, g: u8, n: i64) {
    let g = g as usize;
    LIVE[t][g].fetch_add(n, Relaxed);
    CUM[t][g].fetch_add(n as u64, Relaxed);
    CNT[t][g].fetch_add(1, Relaxed);
    let tot = TOTAL.fetch_add(n, Relaxed) + n;
    if tot > PEAK.load(Relaxed) {
        PEAK.fetch_max(tot, Relaxed);
        if tot >= SNAP_AT.load(Relaxed) {
            snapshot(tot);
        }
    }
}
#[inline]
fn sub(t: usize, g: u8, n: i64) {
    LIVE[t][g as usize].fetch_sub(n, Relaxed);
    TOTAL.fetch_sub(n, Relaxed);
}

pub struct ProbeAlloc;

#[inline]
fn hdr(align: usize) -> usize {
    align.max(16)
}

unsafe impl GlobalAlloc for ProbeAlloc {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let h = hdr(l.align());
        let base = System.alloc(Layout::from_size_align_unchecked(l.size() + h, l.align()));
        if base.is_null() {
            return base;
        }
        let p = base.add(h);
        let (t, g) = (tid(), cur_tag());
        (p.sub(16) as *mut [u32; 4]).write([t as u32, g as u32, 0, 0]);
        add(t, g, l.size() as i64);
        p
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        let h = hdr(l.align());
        let base = System.alloc_zeroed(Layout::from_size_align_unchecked(l.size() + h, l.align()));
        if base.is_null() {
            return base;
        }
        let p = base.add(h);
        let (t, g) = (tid(), cur_tag());
        (p.sub(16) as *mut [u32; 4]).write([t as u32, g as u32, 0, 0]);
        add(t, g, l.size() as i64);
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        let h = hdr(l.align());
        let [t, g, _, _] = (p.sub(16) as *const [u32; 4]).read();
        sub(t as usize, g as u8, l.size() as i64);
        System.dealloc(p.sub(h), Layout::from_size_align_unchecked(l.size() + h, l.align()));
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new_size: usize) -> *mut u8 {
        let h = hdr(l.align());
        let [t, g, _, _] = (p.sub(16) as *const [u32; 4]).read();
        let base = System.realloc(
            p.sub(h),
            Layout::from_size_align_unchecked(l.size() + h, l.align()),
            new_size + h,
        );
        if base.is_null() {
            return base;
        }
        sub(t as usize, g as u8, l.size() as i64);
        let np = base.add(h);
        let (t2, g2) = (tid(), cur_tag());
        (np.sub(16) as *mut [u32; 4]).write([t2 as u32, g2 as u32, 0, 0]);
        add(t2, g2, new_size as i64);
        np
    }
}

/// Move a live allocation (from this allocator, `bytes` = its layout size) to tag `g`.
///
/// # Safety
/// `ptr` must be a live allocation of this allocator of exactly `bytes` bytes.
pub unsafe fn retag(ptr: *const u8, bytes: usize, g: u8) {
    if bytes == 0 {
        return;
    }
    let h = ptr.sub(16) as *mut [u32; 4];
    let [t, og, a, b] = h.read();
    if og as u8 == g {
        return;
    }
    LIVE[t as usize][og as usize].fetch_sub(bytes as i64, Relaxed);
    LIVE[t as usize][g as usize].fetch_add(bytes as i64, Relaxed);
    h.write([t, g as u32, a, b]);
}

/// Retag a Vec's buffer.
pub fn retag_vec<T>(v: &Vec<T>, g: u8) {
    let bytes = v.capacity() * std::mem::size_of::<T>();
    if bytes > 0 {
        unsafe { retag(v.as_ptr() as *const u8, bytes, g) };
    }
}

#[global_allocator]
static GLOBAL: ProbeAlloc = ProbeAlloc;

/// Phase mark (call on the main thread): snapshot of the live table.
pub fn phase(p: usize) {
    let p = p.min(NPH - 1);
    PHASE.store(p as u8, Relaxed);
    copy_tab(&PHASE_SNAP[p]);
    let (_, anon) = statm();
    let v = [1, now_us(), TOTAL.load(Relaxed) as u64, anon];
    for (a, b) in PHASE_INFO[p].iter().zip(v) {
        a.store(b, Relaxed);
    }
}

static STARTED: AtomicBool = AtomicBool::new(false);
static STOP: AtomicBool = AtomicBool::new(false);

/// Start the RssAnon sampler (2 ms).
pub fn start() {
    let mut ts = [0i64; 2];
    unsafe { clock_gettime(1, &mut ts) };
    T0_NS.store(ts[0] as u64 * 1_000_000_000 + ts[1] as u64, Relaxed);
    let _ = Instant::now();
    if STARTED.swap(true, Relaxed) {
        return;
    }
    let _ = std::thread::Builder::new().name("probe-rss".into()).spawn(|| {
        let mut best = 0u64;
        while !STOP.load(Relaxed) {
            let (res, anon) = statm();
            if anon > best {
                best = anon;
                copy_tab(&RSS_SNAP);
                let mi = unsafe { mallinfo2() };
                let v = [
                    anon,
                    PHASE.load(Relaxed) as u64,
                    now_us(),
                    TOTAL.load(Relaxed) as u64,
                    res,
                    mi.arena as u64,
                    mi.uordblks as u64,
                    mi.fordblks as u64,
                ];
                for (a, b) in RSS_INFO.iter().zip(v) {
                    a.store(b, Relaxed);
                }
            }
            std::thread::sleep(std::time::Duration::from_millis(2));
        }
    });
}

// ---- object counters (trees, scratches), per thread ----
pub const NOBJ: usize = 12;
pub const OBJ_NAMES: [&str; NOBJ] = [
    "tree_new", "tree_live_max", "tree_nodes_cap_bytes_sum", "tree_nodes_cap_bytes_max",
    "tree_grid_bytes_sum", "tree_grid_bytes_max", "tree_ends_bytes_sum", "tree_other_bytes_max",
    "grid_allocs", "scratch_new", "scratch_live_max", "tree_drop",
];
#[allow(clippy::declare_interior_mutable_const)]
const ROWO: [AtomicU64; NOBJ] = [ZU; NOBJ];
static OBJ: [[AtomicU64; NOBJ]; NT] = [ROWO; NT];
static TREE_LIVE: [AtomicI64; NT] = [ZI; NT];
static SCR_LIVE: [AtomicI64; NT] = [ZI; NT];

pub fn obj_add(i: usize, v: u64) {
    OBJ[tid()][i].fetch_add(v, Relaxed);
}
pub fn obj_max(i: usize, v: u64) {
    OBJ[tid()][i].fetch_max(v, Relaxed);
}
pub fn tree_created() {
    let t = tid();
    OBJ[t][0].fetch_add(1, Relaxed);
    let l = TREE_LIVE[t].fetch_add(1, Relaxed) + 1;
    OBJ[t][1].fetch_max(l as u64, Relaxed);
}
pub fn tree_dropped(nodes: u64, grid: u64, ends: u64, other: u64) {
    let t = tid();
    TREE_LIVE[t].fetch_sub(1, Relaxed);
    let o = &OBJ[t];
    o[11].fetch_add(1, Relaxed);
    o[2].fetch_add(nodes, Relaxed);
    o[3].fetch_max(nodes, Relaxed);
    o[4].fetch_add(grid, Relaxed);
    o[5].fetch_max(grid, Relaxed);
    o[6].fetch_add(ends, Relaxed);
    o[7].fetch_max(other, Relaxed);
}
pub fn scratch_created() {
    let t = tid();
    OBJ[t][9].fetch_add(1, Relaxed);
    let l = SCR_LIVE[t].fetch_add(1, Relaxed) + 1;
    OBJ[t][10].fetch_max(l as u64, Relaxed);
}
pub fn scratch_dropped() {
    SCR_LIVE[tid()].fetch_sub(1, Relaxed);
}

fn write_tab(w: &mut dyn Write, name: &str, tab: &[[AtomicI64; NG]; NT]) {
    let _ = writeln!(w, "## {name} (KiB live by allocating thread x tag; tid 0 = main)");
    let mut sums = [0i64; NG];
    for (t, row) in tab.iter().enumerate() {
        let v: Vec<i64> = row.iter().map(|a| a.load(Relaxed)).collect();
        if v.iter().all(|&x| x.abs() < 1024) {
            continue;
        }
        let tot: i64 = v.iter().sum();
        let mut s = format!("t{t:<3} total {:>8}", tot >> 10);
        for g in 0..NG {
            sums[g] += v[g];
            if v[g].abs() >= 1024 {
                s.push_str(&format!(" {}={}", TAG_NAMES[g], v[g] >> 10));
            }
        }
        let _ = writeln!(w, "{s}");
    }
    let mut s = format!("SUM  total {:>8}", sums.iter().sum::<i64>() >> 10);
    for g in 0..NG {
        if sums[g].abs() >= 1024 {
            s.push_str(&format!(" {}={}", TAG_NAMES[g], sums[g] >> 10));
        }
    }
    let _ = writeln!(w, "{s}");
}

pub fn report() {
    STOP.store(true, Relaxed);
    phase(6);
    let mut out: Box<dyn Write> = match std::env::var("PROBE_MEM_OUT") {
        Ok(p) => match std::fs::File::create(p) {
            Ok(f) => Box::new(f),
            Err(_) => Box::new(std::io::stderr()),
        },
        Err(_) => Box::new(std::io::stderr()),
    };
    let w: &mut dyn Write = &mut *out;
    let si: Vec<u64> = SNAP_INFO.iter().map(|a| a.load(Relaxed)).collect();
    let _ = writeln!(
        w,
        "# live peak {} KiB (snap at {} KiB) phase {} t {:.3}s mallinfo2 arena {} uord {} ford {} hblkhd {} KiB anon {} KiB",
        PEAK.load(Relaxed) >> 10,
        si[0] >> 10,
        PHASE_NAMES[si[1] as usize],
        si[2] as f64 / 1e6,
        si[3] >> 10,
        si[4] >> 10,
        si[5] >> 10,
        si[6] >> 10,
        si[7]
    );
    write_tab(w, "LIVE PEAK", &SNAP);
    let ri: Vec<u64> = RSS_INFO.iter().map(|a| a.load(Relaxed)).collect();
    let _ = writeln!(
        w,
        "# anon peak {} KiB (resident {} KiB) phase {} t {:.3}s live {} KiB mallinfo2 arena {} uord {} ford {} KiB",
        ri[0],
        ri[4],
        PHASE_NAMES[ri[1] as usize],
        ri[2] as f64 / 1e6,
        ri[3] >> 10,
        ri[5] >> 10,
        ri[6] >> 10,
        ri[7] >> 10
    );
    write_tab(w, "ANON PEAK", &RSS_SNAP);
    for p in 0..NPH {
        let pi: Vec<u64> = PHASE_INFO[p].iter().map(|a| a.load(Relaxed)).collect();
        if pi[0] == 0 {
            continue;
        }
        let _ = writeln!(
            w,
            "# phase {} t {:.3}s live {} KiB anon {} KiB",
            PHASE_NAMES[p],
            pi[1] as f64 / 1e6,
            pi[2] >> 10,
            pi[3]
        );
        write_tab(w, &format!("PHASE {}", PHASE_NAMES[p]), &PHASE_SNAP[p]);
    }
    let _ = writeln!(w, "## cumulative MiB allocated / count by tag (all threads; main separately)");
    for g in 0..NG {
        let (mut cm, mut cp, mut nm, mut np) = (0u64, 0u64, 0u64, 0u64);
        for t in 0..NT {
            let c = CUM[t][g].load(Relaxed);
            let n = CNT[t][g].load(Relaxed);
            if t == 0 {
                cm += c;
                nm += n;
            } else {
                cp += c;
                np += n;
            }
        }
        if nm + np > 0 {
            let _ = writeln!(
                w,
                "{:<12} main {:>9} MiB {:>10} allocs | others {:>9} MiB {:>10} allocs",
                TAG_NAMES[g],
                cm >> 20,
                nm,
                cp >> 20,
                np
            );
        }
    }
    let _ = writeln!(w, "## objects per thread");
    for t in 0..NT {
        let v: Vec<u64> = OBJ[t].iter().map(|a| a.load(Relaxed)).collect();
        if v.iter().all(|&x| x == 0) {
            continue;
        }
        let mut s = format!("t{t:<3}");
        for i in 0..NOBJ {
            if v[i] != 0 {
                s.push_str(&format!(" {}={}", OBJ_NAMES[i], v[i]));
            }
        }
        let _ = writeln!(w, "{s}");
    }
    let _ = w.flush();
}
