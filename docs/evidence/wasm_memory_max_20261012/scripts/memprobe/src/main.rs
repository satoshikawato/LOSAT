//! Probe of a threaded wasm32 WASI module near the top of its linear memory.
//! Usage: memprobe <block MiB> <threads> <scratch file> [oom]
//! It allocates blocks until the allocator refuses, tags one word per MiB, then checks,
//! above the highest block reached: bulk fill/copy, threads (stack address, clock,
//! random, mutex/condvar, stdout), and file write/read with high buffers.
use std::alloc::{alloc, dealloc, Layout};
use std::collections::HashMap;
use std::io::{Read, Write};
use std::sync::{Arc, Condvar, Mutex};
use std::time::Instant;

const MIB: usize = 1 << 20;

fn pages() -> usize {
    #[cfg(target_arch = "wasm32")]
    {
        core::arch::wasm32::memory_size::<0>()
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        0
    }
}

fn tag(block: usize, mib: usize) -> u64 {
    ((block as u64) << 32 | mib as u64) ^ 0x5a5a_a5a5_c3c3_3c3c
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let block = args.get(1).and_then(|s| s.parse::<usize>().ok()).unwrap_or(256) * MIB;
    let threads = args.get(2).and_then(|s| s.parse::<usize>().ok()).unwrap_or(4);
    let scratch = args.get(3).cloned().unwrap_or_default();
    let oom = args.get(4).map(|s| s == "oom").unwrap_or(false);
    println!("start pages={} block_mib={}", pages(), block / MIB);

    let layout = Layout::from_size_align(block, 65536).unwrap();
    // The bookkeeping is reserved first: once memory is full it cannot grow.
    let mut blocks: Vec<*mut u8> = Vec::with_capacity(8192);
    let mut small: Vec<(*mut u8, Layout)> = Vec::with_capacity(3_000_000);
    loop {
        let p = unsafe { alloc(layout) };
        if p.is_null() {
            break;
        }
        for mib in 0..block / MIB {
            unsafe { (p.add(mib * MIB) as *mut u64).write(tag(blocks.len(), mib)) };
        }
        unsafe { (p.add(block - 8) as *mut u64).write(tag(blocks.len(), 0xfff_fff)) };
        blocks.push(p);
    }
    let lowest = blocks.iter().map(|p| *p as usize).min().unwrap_or(0);
    let highest = blocks.iter().map(|p| *p as usize).max().unwrap_or(0);
    println!(
        "blocks={} total_mib={} lowest=0x{:08x} highest_start=0x{:08x} highest_end=0x{:08x} pages={}",
        blocks.len(),
        blocks.len() * (block / MIB),
        lowest,
        highest,
        highest.wrapping_add(block),
        pages()
    );
    let mut bad = 0usize;
    for (i, p) in blocks.iter().enumerate() {
        for mib in 0..block / MIB {
            if unsafe { (p.add(mib * MIB) as *const u64).read() } != tag(i, mib) {
                bad += 1;
            }
        }
        if unsafe { (p.add(block - 8) as *const u64).read() } != tag(i, 0xfff_fff) {
            bad += 1;
        }
    }
    println!("tag_mismatches={}", bad);

    // Bulk operations at the top block (memory.fill / memory.copy with +bulk-memory).
    let top = blocks.iter().copied().max_by_key(|p| *p as usize).unwrap();
    let n = (3 * MIB + 12345).min(block / 4);
    unsafe {
        std::ptr::write_bytes(top.add(block - n), 0xab, n);
        std::ptr::copy_nonoverlapping(top.add(block - n), top.add(block - 2 * n - 7), n);
        let a = std::slice::from_raw_parts(top.add(block - n), n);
        let b = std::slice::from_raw_parts(top.add(block - 2 * n - 7), n);
        println!("bulk_ok={} last_byte_addr=0x{:08x}", a == b && a.iter().all(|x| *x == 0xab), top as usize + block - 1);
    }

    // Fill what is left with small chunks, down to the last bytes below the maximum.
    for size in [60 * 1024usize, 1024, 16] {
        let l = Layout::from_size_align(size, 8).unwrap();
        loop {
            if small.len() == small.capacity() {
                break;
            }
            let p = unsafe { alloc(l) };
            if p.is_null() {
                break;
            }
            unsafe { (p as *mut u64).write(tag(small.len(), size)) };
            small.push((p, l));
            if small.len() % 100_000 == 0 || (std::env::var_os("PROBE_TRACE").is_some() && small.len() % 50 == 0) {
                eprintln!("fill size={} chunks={} last=0x{:08x} pages={}", size, small.len(), p as usize, pages());
            }
        }
        eprintln!("fill size={} done chunks={} pages={}", size, small.len(), pages());
    }
    let small_top = small.iter().map(|(p, l)| *p as usize as u64 + l.size() as u64).max().unwrap_or(0);
    let small_bad = small
        .iter()
        .enumerate()
        .filter(|(i, (p, l))| unsafe { (*p as *const u64).read() } != tag(*i, l.size()))
        .count();
    println!("small_chunks={} small_top=0x{:09x} small_mismatches={} pages={}", small.len(), small_top, small_bad, pages());
    for (p, l) in small.drain(..) {
        unsafe { dealloc(p, l) };
    }

    // Free the highest blocks (at least 96 MiB) and one low block so that threads and
    // buffers find memory, high first.
    blocks.sort_by_key(|p| *p as usize);
    let low = blocks.remove(0);
    let mut freed = 0usize;
    while freed < 96 * MIB {
        let p = blocks.pop().expect("a block to free");
        unsafe { dealloc(p, layout) };
        freed += block;
    }
    let high_buf_layout = Layout::from_size_align(16 * MIB, 65536).unwrap();
    let hold = unsafe { alloc(high_buf_layout) };
    assert!(!hold.is_null(), "no memory for the high buffer");
    println!("hold=0x{:08x}", hold as usize);

    // Threads: stacks and their WASI calls.
    let shared = Arc::new((Mutex::new(0usize), Condvar::new()));
    let start = Instant::now();
    let mut handles = Vec::new();
    for t in 0..threads {
        let shared = Arc::clone(&shared);
        handles.push(
            std::thread::Builder::new()
                .stack_size(4 * MIB)
                .spawn(move || {
                    let local = 0u8;
                    let stack = &local as *const u8 as usize;
                    let now = Instant::now();
                    let mut map: HashMap<u64, u64> = HashMap::new();
                    for i in 0..1000u64 {
                        map.insert(i, i * i);
                    }
                    let v: Vec<u8> = vec![t as u8; 3 * MIB];
                    let sum: u64 = v.iter().map(|x| *x as u64).sum::<u64>() + map[&999];
                    let (lock, cv) = &*shared;
                    let mut g = lock.lock().unwrap();
                    *g += 1;
                    cv.notify_all();
                    drop(g);
                    println!(
                        "thread {} stack=0x{:08x} vec=0x{:08x} sum={} elapsed_ns>0={}",
                        t,
                        stack,
                        v.as_ptr() as usize,
                        sum,
                        now.elapsed().as_nanos() > 0
                    );
                    stack
                })
                .expect("spawn"),
        );
    }
    {
        let (lock, cv) = &*shared;
        let mut g = lock.lock().unwrap();
        while *g < threads {
            g = cv.wait(g).unwrap();
        }
    }
    let stacks: Vec<usize> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    println!(
        "threads_joined={} stacks_above_2g={} elapsed_ok={}",
        stacks.len(),
        stacks.iter().filter(|s| **s >= 0x8000_0000).count(),
        start.elapsed().as_nanos() > 0
    );

    // File write and read with high buffers.
    if !scratch.is_empty() {
        let size = 5 * MIB + 321;
        let a = unsafe { std::slice::from_raw_parts_mut(hold, size) };
        for (i, x) in a.iter_mut().enumerate() {
            *x = (i * 31 + 7) as u8;
        }
        std::fs::File::create(&scratch).and_then(|mut f| f.write_all(a)).expect("write scratch");
        let b = unsafe { std::slice::from_raw_parts_mut(hold.add(8 * MIB), size) };
        let mut f = std::fs::File::open(&scratch).expect("open scratch");
        f.read_exact(b).expect("read scratch");
        println!("file_ok={} buf=0x{:08x}", a == b, hold as usize);
        let _ = std::fs::remove_file(&scratch);
    }
    unsafe { dealloc(hold, high_buf_layout) };
    // The blocks that stayed allocated through the fill, the threads and the file I/O
    // must still hold their tags (block indices are the addresses' ranks at tagging time,
    // so compare the per-MiB part only).
    let mut late_bad = 0usize;
    for p in blocks.iter().chain(std::iter::once(&low)) {
        for mib in 0..block / MIB {
            let word = unsafe { (p.add(mib * MIB) as *const u64).read() } ^ 0x5a5a_a5a5_c3c3_3c3c;
            if (word & 0xffff_ffff) as usize != mib {
                late_bad += 1;
            }
        }
    }
    println!("late_tag_mismatches={}", late_bad);
    unsafe { dealloc(low, layout) };
    for p in blocks {
        unsafe { dealloc(p, layout) };
    }
    println!("end pages={}", pages());
    if oom {
        // The infallible path: Rust prints "memory allocation of N bytes failed" and aborts.
        let mut keep: Vec<Vec<u8>> = Vec::new();
        loop {
            keep.push(Vec::with_capacity(block));
        }
    }
}
