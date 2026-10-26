//! The global allocator of the threaded Wasm builds (wasm32 with the `wasm-threads` feature: the
//! threaded command and reactor): caches of small free blocks, used by one thread at a time, in
//! front of the system allocator (wasi-libc's malloc).
//!
//! No NCBI counterpart: this decides only where a freed block waits until it is used again; no
//! value that LOSAT computes changes. NCBI allocates with the C library's malloc, which on Linux
//! (glibc) keeps a cache of small free blocks per thread (tcache) and gives threads arenas of
//! their own, so threads seldom wait for each other's allocations. wasi-libc's malloc (dlmalloc)
//! has one heap behind one lock, which a waiting thread spins for, calling `sched_yield` every
//! 64 tries. LOSAT's searches allocate and free often (BLASTN megablast on one thread: 7.6
//! million allocations in 8.8 s, all of at most 64 KiB), so at 16 threads most of the time went
//! to that lock: in the threaded command's megablast of 100 queries x 16S sequences, 75 % of the
//! CPU samples had dlmalloc, dlfree or realloc on the stack (BLASTN: 16 %), and 16 threads ran
//! 2.5 times faster than one (native: 14.6).
//!
//! A cache keeps, per size class of at most `LARGEST` bytes and per alignment (8, or 16), a list
//! of free blocks of the class's size, each allocated from the system allocator with the class's
//! layout, up to `CAPS` blocks per class and `CACHE_BYTES` in all; a block beyond those goes back
//! to the system allocator. A thread locks a cache for each allocation or free (a lock that is
//! taken goes to the next cache instead), and starts from the cache that it used last. A thread
//! that keeps finding its cache last used by another thread moves on to the next one, so the
//! threads that run at the same time settle on caches of their own, the lowest ones free (a new
//! thread starts at the first). A thread keeps no cache when it ends (a global allocator may
//! not use thread-local values with drops, and LOSAT's code calls no C function): the next
//! thread uses that cache's blocks. So the blocks held are bounded by `CACHES_LEN` caches, and in
//! practice by the threads that run at once. A larger or more aligned request goes to the system
//! allocator as before.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::{Cell, UnsafeCell};
use std::ptr;
use std::sync::atomic::{AtomicBool, Ordering};

#[cfg(all(target_arch = "wasm32", feature = "wasm-threads"))]
#[global_allocator]
static ALLOCATOR: ThreadCacheAlloc = ThreadCacheAlloc;

/// The caches (a power of two).
const CACHES_LEN: usize = 64;
/// The largest cached size.
const LARGEST: usize = 64 << 10;
/// Size classes: 16 to 256 bytes in steps of 16, then four per doubling up to `LARGEST`.
const CLASSES: usize = 16 + 4 * (LARGEST.trailing_zeros() as usize - 8);
/// The alignments of the cached blocks: 8 (the system allocator's malloc on wasm32) and 16.
const ALIGNS: [usize; 2] = [8, 16];
/// The most bytes that one cache holds.
const CACHE_BYTES: usize = 1 << 20;
/// The most bytes of one class and alignment that one cache holds (at least 4 blocks).
const CLASS_BYTES: usize = 64 << 10;
/// A thread moves to the next cache when this many of its operations in a row found the cache
/// last used by another thread.
const MOVE_AFTER: u32 = 4;

const fn class_size(class: usize) -> usize {
    if class < 16 {
        (class + 1) << 4
    } else {
        let lg = (class - 16) / 4 + 8;
        (1 << lg) + (((class - 16) % 4 + 1) << (lg - 2))
    }
}

const SIZES: [usize; CLASSES] = {
    let mut sizes = [0; CLASSES];
    let mut class = 0;
    while class < CLASSES {
        sizes[class] = class_size(class);
        class += 1;
    }
    sizes
};

const CAPS: [u32; CLASSES] = {
    let mut caps = [0; CLASSES];
    let mut class = 0;
    while class < CLASSES {
        let blocks = CLASS_BYTES / class_size(class);
        caps[class] = if blocks < 4 { 4 } else { blocks as u32 };
        class += 1;
    }
    caps
};

/// The class of the smallest size class that holds `size` bytes (`size <= LARGEST`).
#[inline(always)]
fn class_of(size: usize) -> usize {
    if size <= 256 {
        size.saturating_sub(1) >> 4
    } else {
        let s = size - 1;
        let lg = (usize::BITS - 1 - s.leading_zeros()) as usize;
        16 + (lg - 8) * 4 + ((s >> (lg - 2)) & 3)
    }
}

/// The class and alignment (index into `ALIGNS`) of a request that the caches serve, if any.
#[inline(always)]
fn small(layout: Layout) -> Option<(usize, usize)> {
    (layout.size() <= LARGEST && layout.align() <= ALIGNS[1]).then(|| {
        (
            class_of(layout.size()),
            usize::from(layout.align() > ALIGNS[0]),
        )
    })
}

/// The layout with which the blocks of a class and alignment are allocated and freed.
#[inline(always)]
fn class_layout(class: usize, align: usize) -> Layout {
    // SAFETY: every class size is a non-zero multiple of 16, and both alignments are powers of
    // two; the size is far below isize::MAX.
    unsafe { Layout::from_size_align_unchecked(SIZES[class], ALIGNS[align]) }
}

struct Lists {
    /// The first free block of each alignment and class; a free block's first word holds the
    /// next one.
    head: [[*mut u8; CLASSES]; 2],
    count: [[u32; CLASSES]; 2],
    bytes: usize,
    /// The thread that used the cache last (`token`).
    last: usize,
}

impl Lists {
    /// A free block of the class and alignment, or null.
    #[inline(always)]
    fn pop(&mut self, class: usize, align: usize) -> *mut u8 {
        let block = self.head[align][class];
        if !block.is_null() {
            // SAFETY: a listed block is a free block of at least 16 bytes, aligned to at least 8,
            // whose first word push wrote.
            self.head[align][class] = unsafe { *(block as *mut *mut u8) };
            self.count[align][class] -= 1;
            self.bytes -= SIZES[class];
        }
        block
    }

    /// Keeps `block`, unless the cache is full.
    ///
    /// # Safety
    /// `block` is a free block allocated with `class_layout(class, align)`.
    #[inline(always)]
    unsafe fn push(&mut self, class: usize, align: usize, block: *mut u8) -> bool {
        if self.count[align][class] >= CAPS[class] || self.bytes + SIZES[class] > CACHE_BYTES {
            return false;
        }
        // SAFETY: the block is free, at least 16 bytes and aligned to at least 8 (the caller's
        // contract).
        unsafe { *(block as *mut *mut u8) = self.head[align][class] };
        self.head[align][class] = block;
        self.count[align][class] += 1;
        self.bytes += SIZES[class];
        true
    }
}

#[repr(align(64))]
struct Cache {
    locked: AtomicBool,
    lists: UnsafeCell<Lists>,
}

// SAFETY: `lists` is used only by the thread that set `locked` (with Acquire, cleared with
// Release, so the next one sees the lists as the last one left them).
unsafe impl Sync for Cache {}

impl Cache {
    const fn new() -> Self {
        Self {
            locked: AtomicBool::new(false),
            lists: UnsafeCell::new(Lists {
                head: [[ptr::null_mut(); CLASSES]; 2],
                count: [[0; CLASSES]; 2],
                bytes: 0,
                last: 0,
            }),
        }
    }
}

static CACHES: [Cache; CACHES_LEN] = [const { Cache::new() }; CACHES_LEN];

thread_local! {
    /// The cache that the thread tries first.
    static HINT: Cell<usize> = const { Cell::new(0) };
    /// How many of the thread's operations in a row found the cache last used by another thread.
    static STRIKES: Cell<u32> = const { Cell::new(0) };
}

/// A number that no other running thread has: the address of its `HINT`.
#[inline(always)]
fn token() -> usize {
    HINT.with(|hint| hint as *const Cell<usize> as usize)
}

/// Runs `f` on the lists of a cache that the calling thread locked: the one it used last, else
/// the next one free; `None` if every cache was locked.
#[inline(always)]
fn with_cache<R>(f: impl FnOnce(&mut Lists) -> R) -> Option<R> {
    let hint = HINT.get();
    for probe in 0..CACHES_LEN {
        let index = (hint + probe) % CACHES_LEN;
        let cache = &CACHES[index];
        if cache.locked.swap(true, Ordering::Acquire) {
            continue;
        }
        // SAFETY: this thread holds the cache's lock.
        let lists = unsafe { &mut *cache.lists.get() };
        let mut next = index;
        let me = token();
        if lists.last == me {
            STRIKES.set(0);
        } else {
            lists.last = me;
            let strikes = STRIKES.get() + 1;
            if strikes >= MOVE_AFTER {
                next = (index + 1) % CACHES_LEN;
                STRIKES.set(0);
            } else {
                STRIKES.set(strikes);
            }
        }
        let result = f(lists);
        cache.locked.store(false, Ordering::Release);
        if next != hint {
            HINT.set(next);
        }
        return Some(result);
    }
    None
}

struct ThreadCacheAlloc;

// SAFETY: every block handed out was allocated by System with a layout of at least the request's
// size and alignment (a class block: the class's size and an alignment that covers the
// request's), and every block handed back goes to System with the layout it was allocated with
// (a class block: its class's), or to a cache that hands it out only once again.
unsafe impl GlobalAlloc for ThreadCacheAlloc {
    #[inline]
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if let Some((class, align)) = small(layout) {
            if let Some(block) = with_cache(|lists| lists.pop(class, align)) {
                if !block.is_null() {
                    return block;
                }
            }
            // SAFETY: a class layout has a non-zero size.
            return unsafe { System.alloc(class_layout(class, align)) };
        }
        // SAFETY: the caller's layout, as given.
        unsafe { System.alloc(layout) }
    }

    #[inline]
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        if let Some((class, align)) = small(layout) {
            if let Some(block) = with_cache(|lists| lists.pop(class, align)) {
                if !block.is_null() {
                    // SAFETY: the block holds at least layout.size() bytes.
                    unsafe { ptr::write_bytes(block, 0, layout.size()) };
                    return block;
                }
            }
            // SAFETY: a class layout has a non-zero size.
            return unsafe { System.alloc_zeroed(class_layout(class, align)) };
        }
        // SAFETY: the caller's layout, as given.
        unsafe { System.alloc_zeroed(layout) }
    }

    #[inline]
    unsafe fn dealloc(&self, block: *mut u8, layout: Layout) {
        if let Some((class, align)) = small(layout) {
            // SAFETY: a block of this layout was allocated with the class's layout (alloc and
            // realloc here).
            if with_cache(|lists| unsafe { lists.push(class, align, block) }) == Some(true) {
                return;
            }
            // SAFETY: the block was allocated by System with the class's layout.
            return unsafe { System.dealloc(block, class_layout(class, align)) };
        }
        // SAFETY: the block was allocated by System with this layout.
        unsafe { System.dealloc(block, layout) }
    }

    #[inline]
    unsafe fn realloc(&self, block: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        // SAFETY: the caller guarantees that new_size, rounded up to the alignment, does not
        // overflow isize.
        let new_layout = unsafe { Layout::from_size_align_unchecked(new_size, layout.align()) };
        match (small(layout), small(new_layout)) {
            // SAFETY: the block was allocated by System with this layout.
            (None, None) => unsafe { System.realloc(block, layout, new_size) },
            // The block already holds the class's size.
            (Some(old), Some(new)) if old == new => block,
            _ => {
                // SAFETY: new_layout has a non-zero size (the caller's contract).
                let moved = unsafe { self.alloc(new_layout) };
                if !moved.is_null() {
                    // SAFETY: both blocks hold at least the bytes copied, and they are distinct.
                    unsafe {
                        ptr::copy_nonoverlapping(block, moved, layout.size().min(new_size));
                        self.dealloc(block, layout);
                    }
                }
                moved
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;

    #[test]
    fn classes_hold_their_sizes() {
        assert_eq!(SIZES[CLASSES - 1], LARGEST);
        for size in 1..=LARGEST {
            let class = class_of(size);
            assert!(SIZES[class] >= size, "size {size} class {class}");
            assert!(
                class == 0 || SIZES[class - 1] < size,
                "size {size} class {class}"
            );
        }
        assert!(SIZES.iter().all(|size| size % 16 == 0));
        assert!(CAPS.iter().all(|&caps| caps >= 4));
    }

    struct Block {
        at: *mut u8,
        layout: Layout,
        fill: u8,
    }

    // SAFETY (test): a block is used by one thread at a time; it moves to another only to be freed.
    unsafe impl Send for Block {}

    fn check(block: &Block) {
        for i in 0..block.layout.size() {
            // SAFETY: the block holds layout.size() bytes.
            assert_eq!(
                unsafe { *block.at.add(i) },
                block.fill,
                "byte {i} of {:?}",
                block.layout
            );
        }
    }

    fn new_block(size: usize, align: usize, fill: u8, zeroed: bool) -> Block {
        let layout = Layout::from_size_align(size, align).unwrap();
        // SAFETY: a non-zero size.
        let at = unsafe {
            if zeroed {
                ThreadCacheAlloc.alloc_zeroed(layout)
            } else {
                ThreadCacheAlloc.alloc(layout)
            }
        };
        assert!(!at.is_null());
        assert_eq!(at as usize % align, 0);
        if zeroed {
            check(&Block {
                at,
                layout,
                fill: 0,
            });
        }
        // SAFETY: the block holds `size` bytes.
        unsafe { ptr::write_bytes(at, fill, size) };
        Block { at, layout, fill }
    }

    fn release(block: Block) {
        check(&block);
        // SAFETY: the block came from ThreadCacheAlloc with this layout.
        unsafe { ThreadCacheAlloc.dealloc(block.at, block.layout) };
    }

    // Every cache is unlocked, within its bounds, and its lists hold as many blocks as counted,
    // each aligned to its list's alignment.
    fn check_caches() {
        for cache in &CACHES {
            assert!(!cache.locked.swap(true, Ordering::Acquire));
            // SAFETY: the lock is held.
            let lists = unsafe { &*cache.lists.get() };
            let mut bytes = 0;
            for align in 0..2 {
                for class in 0..CLASSES {
                    let (mut block, mut blocks) = (lists.head[align][class], 0);
                    while !block.is_null() {
                        assert_eq!(block as usize % ALIGNS[align], 0);
                        blocks += 1;
                        // SAFETY: a listed block holds the next one in its first word.
                        block = unsafe { *(block as *mut *mut u8) };
                    }
                    assert_eq!(blocks, lists.count[align][class]);
                    assert!(blocks <= CAPS[class]);
                    bytes += blocks as usize * SIZES[class];
                }
            }
            assert_eq!(bytes, lists.bytes);
            assert!(bytes <= CACHE_BYTES);
            cache.locked.store(false, Ordering::Release);
        }
    }

    // Blocks of every class and beyond, zeroed or not, more aligned, grown and shrunk across
    // classes, freed on their own thread or on another: every byte stays as written, and the
    // caches stay consistent.
    #[test]
    fn blocks_keep_their_bytes_across_threads() {
        let threads = 6;
        let (senders, receivers): (Vec<_>, Vec<_>) =
            (0..threads).map(|_| mpsc::channel::<Block>()).unzip();
        std::thread::scope(|scope| {
            for (t, inbox) in receivers.into_iter().enumerate() {
                let next = senders[(t + 1) % threads].clone();
                scope.spawn(move || {
                    let mut seed = 0x9e37_79b9_u32.wrapping_mul(t as u32 + 1);
                    let mut rand = move || {
                        seed ^= seed << 13;
                        seed ^= seed >> 17;
                        seed ^= seed << 5;
                        seed as usize
                    };
                    let mut live: Vec<Block> = Vec::new();
                    for step in 0..4000 {
                        let size = match rand() % 8 {
                            0 => 1 + rand() % 16,
                            1..=4 => 1 + rand() % 4096,
                            5 | 6 => 1 + rand() % LARGEST,
                            _ => LARGEST + 1 + rand() % 50_000,
                        };
                        let align = [1, 8, 16, 64][rand() % 4];
                        live.push(new_block(size, align, step as u8, rand() % 3 == 0));
                        if rand() % 4 == 0 {
                            let at = rand() % live.len();
                            let old = &mut live[at];
                            check(old);
                            let new_size = 1 + rand() % (2 * LARGEST);
                            // SAFETY: the block came from ThreadCacheAlloc with old.layout.
                            let moved =
                                unsafe { ThreadCacheAlloc.realloc(old.at, old.layout, new_size) };
                            assert!(!moved.is_null());
                            let kept = Layout::from_size_align(old.layout.size().min(new_size), 1);
                            check(&Block {
                                at: moved,
                                layout: kept.unwrap(),
                                fill: old.fill,
                            });
                            old.at = moved;
                            old.layout =
                                Layout::from_size_align(new_size, old.layout.align()).unwrap();
                            // SAFETY: the block holds new_size bytes.
                            unsafe { ptr::write_bytes(moved, old.fill, new_size) };
                        }
                        while live.len() > 64 {
                            let block = live.swap_remove(rand() % live.len());
                            if rand() % 2 == 0 {
                                release(block);
                            } else {
                                next.send(block).unwrap();
                            }
                        }
                        while let Ok(block) = inbox.try_recv() {
                            release(block);
                        }
                    }
                    live.into_iter().for_each(release);
                    drop(next);
                    for block in inbox.iter() {
                        release(block);
                    }
                });
            }
            drop(senders);
        });
        check_caches();
    }

    // Threads that run at the same time and allocate in turn end on caches of their own.
    #[test]
    fn threads_running_together_spread_over_caches() {
        let threads = 4;
        let barrier = std::sync::Barrier::new(threads);
        let hints: Vec<usize> = std::thread::scope(|scope| {
            let workers: Vec<_> = (0..threads)
                .map(|_| {
                    scope.spawn(|| {
                        let layout = Layout::from_size_align(48, 8).unwrap();
                        for _ in 0..2000 {
                            barrier.wait();
                            // SAFETY: a non-zero size; the block goes back with its layout.
                            unsafe {
                                ThreadCacheAlloc.dealloc(ThreadCacheAlloc.alloc(layout), layout)
                            };
                        }
                        HINT.get()
                    })
                })
                .collect();
            workers
                .into_iter()
                .map(|worker| worker.join().unwrap())
                .collect()
        });
        let mut distinct = hints.clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(distinct.len(), threads, "{hints:?}");
    }
}
