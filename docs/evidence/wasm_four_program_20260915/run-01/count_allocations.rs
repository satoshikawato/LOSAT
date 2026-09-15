// NCBI reference: c++/src/algo/blast/core/blast_extend.c:219-224,337-341
// calloc(1, sizeof(BlastInitHitList));
// malloc(MIN_INIT_HITLIST_SIZE * sizeof(BlastInitHSP));
// realloc(match_array, num_avail * sizeof(BlastInitHSP));
// Diagnostic-only native allocator counts, never an adoption timing build.
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
pub struct Counting;
static ALLOC: AtomicUsize = AtomicUsize::new(0);
static REALLOC: AtomicUsize = AtomicUsize::new(0);
static DEALLOC: AtomicUsize = AtomicUsize::new(0);
static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
fn add(n: usize) { let current = LIVE.fetch_add(n, Relaxed) + n; PEAK.fetch_max(current, Relaxed); }
// SAFETY: all memory operations delegate unchanged layouts/pointers to System.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOC.fetch_add(1, Relaxed);
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() { add(layout.size()); }
        pointer
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        ALLOC.fetch_add(1, Relaxed);
        let pointer = unsafe { System.alloc_zeroed(layout) };
        if !pointer.is_null() { add(layout.size()); }
        pointer
    }
    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        DEALLOC.fetch_add(1, Relaxed); LIVE.fetch_sub(layout.size(), Relaxed);
        unsafe { System.dealloc(pointer, layout) };
    }
    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        REALLOC.fetch_add(1, Relaxed);
        let next = unsafe { System.realloc(pointer, layout, size) };
        if !next.is_null() { LIVE.fetch_sub(layout.size(), Relaxed); add(size); }
        next
    }
}
#[global_allocator]
static ALLOCATOR: Counting = Counting;
pub fn report() {
    eprintln!("[ALLOCATOR] alloc={} realloc={} dealloc={} live_requested_bytes={} peak_requested_bytes={}",
        ALLOC.load(Relaxed), REALLOC.load(Relaxed), DEALLOC.load(Relaxed), LIVE.load(Relaxed), PEAK.load(Relaxed));
}
