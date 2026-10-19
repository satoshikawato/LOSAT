//! EXPERIMENT (R10_PROBE; line R's R10 probe of plan v4, `v4-R/DESIGN_R10.md`; not for merge): the
//! serial megablast table fill of `build_mb_lookup` (the LOSAT_X_MBDELAY loop of `x_fill_mb_queued`)
//! measured in its own window. No NCBI counterpart and no change of any value: every arm applies
//! the same update to the same words in the same order (`x_mb_apply`).
//!
//! `R10_PROBE` is a comma-separated list:
//! - `t`: only the line below (the fill's time and RSS).
//! - `ctl`: perf stat control around the fill: writes `enable` / `disable` to the FIFO named by
//!   `R10_CTL` and waits for perf's `ack` on `R10_ACK` (`perf stat -D -1 --control fifo:CTL,ACK`).
//! - `d<N>[h2|hw|bf]`: the fill with the prefetch `N` words ahead (N = 0, 1, 2, 4, 8, 16, 32, 64, 128;
//!   0 = apply each word at once, no prefetch); `h2` = PREFETCHT2, `hw` = PREFETCHW, else
//!   PREFETCHT0 (`d16` is the shape of the default LOSAT_X_MBDELAY loop); `bf` = PREFETCHT0 and
//!   the update without the branch on the old head (`x_apply_bf`: the PV word and the helper count
//!   are always written, with 0 added where the branch would skip them; same values).
//! - `dump=<path>`: before the fill, write each word's (cell, q_off_1) as two little-endian u32,
//!   in the fill's order.
//!
//! With any token, one line on stderr at the end of the fill:
//! `[R10] fill_s=<s> words=<n> rss0_kb= rss1_kb= hwm1_kb= arm=<arm>`.

use super::{
    build_unmasked_ranges, db_word_count_exceeds, x_mb_apply, MaskedInterval, XMbTables,
    BLAST2NA_MASK, PV_ARRAY_MASK,
};
use std::io::{Read, Write};
use std::sync::OnceLock;
use std::time::Instant;

#[derive(Default)]
struct Probe {
    any: bool,
    ctl: bool,
    depth: Option<(usize, u8)>,
    dump: Option<String>,
    arm: String,
}

fn probe() -> &'static Probe {
    static P: OnceLock<Probe> = OnceLock::new();
    P.get_or_init(|| {
        let mut p = Probe::default();
        let v = std::env::var("R10_PROBE").unwrap_or_default();
        for tok in v.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            p.any = true;
            let (k, val) = tok.split_once('=').unwrap_or((tok, ""));
            if k == "t" {
                // the time and RSS line only
            } else if k == "ctl" {
                p.ctl = true;
            } else if k == "dump" {
                p.dump = Some(val.to_string());
            } else if let Some(rest) = k.strip_prefix('d') {
                let (n, h) = if let Some(n) = rest.strip_suffix("bf") {
                    (n, 4u8)
                } else if let Some(n) = rest.strip_suffix("h2") {
                    (n, 2u8)
                } else if let Some(n) = rest.strip_suffix("hw") {
                    (n, 3u8)
                } else {
                    (rest, 0u8)
                };
                match n.parse::<usize>() {
                    Ok(d) if [0, 1, 2, 4, 8, 16, 32, 64, 128].contains(&d) => {
                        p.depth = Some((d, h));
                        p.arm = tok.to_string();
                    }
                    _ => eprintln!("[R10] unknown R10_PROBE token {tok}"),
                }
            } else {
                eprintln!("[R10] unknown R10_PROBE token {tok}");
            }
        }
        if p.arm.is_empty() {
            p.arm = "default".to_string();
        }
        p
    })
}

fn status_kb(key: &str) -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with(key))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|x| x.parse().ok())
        })
        .unwrap_or(0)
}

struct Ctl {
    ctl: std::fs::File,
    ack: std::fs::File,
}

fn ctl() -> Option<&'static std::sync::Mutex<Ctl>> {
    static C: OnceLock<Option<std::sync::Mutex<Ctl>>> = OnceLock::new();
    C.get_or_init(|| {
        let c = std::env::var("R10_CTL").ok()?;
        let a = std::env::var("R10_ACK").ok()?;
        let ctl = std::fs::OpenOptions::new().write(true).open(c).ok()?;
        let ack = std::fs::OpenOptions::new().read(true).open(a).ok()?;
        Some(std::sync::Mutex::new(Ctl { ctl, ack }))
    })
    .as_ref()
}

fn send(cmd: &str) {
    let Some(c) = ctl() else {
        eprintln!("[R10] no control fifo (R10_CTL / R10_ACK)");
        return;
    };
    let mut c = c.lock().unwrap();
    let _ = c.ctl.write_all(cmd.as_bytes());
    let _ = c.ctl.flush();
    let mut buf = [0u8; 8];
    let mut got = Vec::new();
    while !got.windows(3).any(|w| w == b"ack") {
        match c.ack.read(&mut buf) {
            Ok(0) | Err(_) => break,
            Ok(n) => got.extend_from_slice(&buf[..n]),
        }
    }
}

pub(super) struct Window {
    t0: Instant,
    rss0: u64,
}

/// Called by `build_mb_lookup` right before the fill (after `x_clock.fill()`).
#[inline(never)]
pub(super) fn begin() -> Option<Window> {
    let p = probe();
    if !p.any {
        return None;
    }
    let rss0 = status_kb("VmRSS:");
    if p.ctl {
        send("enable\n");
    }
    Some(Window {
        t0: Instant::now(),
        rss0,
    })
}

/// Called by `build_mb_lookup` right after the fill (before `x_clock.chain()`).
#[inline(never)]
pub(super) fn end(w: Option<Window>, words: usize) {
    let Some(w) = w else { return };
    let s = w.t0.elapsed().as_secs_f64();
    let p = probe();
    if p.ctl {
        send("disable\n");
    }
    eprintln!(
        "[R10] fill_s={s:.6} words={words} rss0_kb={} rss1_kb={} hwm1_kb={} arm={}",
        w.rss0,
        status_kb("VmRSS:"),
        status_kb("VmHWM:"),
        p.arm
    );
}

/// The probe's fill shape, if one was asked for.
pub(super) fn depth() -> Option<(usize, u8)> {
    probe().depth
}

/// The words of the fill, as `x_fill_mb_queued` reads them.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn for_each_word(
    queries_blastna: &[&[u8]],
    query_offsets: &[i32],
    query_masks: &[Vec<MaskedInterval>],
    word_length: usize,
    lut_word_length: usize,
    kmer_mask: u64,
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
    t: &mut XMbTables<'_>,
    mut f: impl FnMut(&mut XMbTables<'_>, usize, u32),
) {
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna;
        if seq.len() < lut_word_length {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let ranges = build_unmasked_ranges(seq.len(), masks);
        let query_offset = query_offsets[q_idx].max(0) as usize;

        for (range_start, range_end) in ranges {
            let range_len = range_end.saturating_sub(range_start);
            if word_length > range_len {
                continue;
            }

            let mut current_kmer: u64 = 0;
            let mut valid_bases = 0usize;

            for pos in range_start..range_end {
                let base = seq[pos];
                if (base & BLAST2NA_MASK) != 0 {
                    current_kmer = 0;
                    valid_bases = 0;
                    *t.ambiguous_skipped += 1;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | base as u64) & kmer_mask;
                valid_bases += 1;
                if valid_bases < lut_word_length {
                    continue;
                }

                *t.total_positions += 1;

                if let Some(counts_filter) = db_word_counts {
                    if db_word_count_exceeds(counts_filter, current_kmer, max_db_word_count) {
                        continue;
                    }
                }

                let bucket = current_kmer as usize;
                let q_off_1 = (query_offset + (pos + 1 - lut_word_length) + 1) as u32;
                f(t, bucket, q_off_1);
            }
        }
    }
}

#[inline(always)]
fn prefetch<const H: u8>(p: *const u32) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: `p` points into the table (masked k-mer); a prefetch does not fault.
    unsafe {
        match H {
            2 => core::arch::x86_64::_mm_prefetch(p as *const i8, core::arch::x86_64::_MM_HINT_T2),
            3 => {
                core::arch::asm!("prefetchw [{0}]", in(reg) p, options(nostack, preserves_flags, readonly))
            }
            _ => core::arch::x86_64::_mm_prefetch(p as *const i8, core::arch::x86_64::_MM_HINT_T0),
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
}

/// `x_mb_apply` without the branch on the old head: the PV word is ORed with 0 and the helper count
/// gets a saturating 0 where the branch would skip them, so every value is the one `x_mb_apply`
/// writes.
#[inline(always)]
fn x_apply_bf(t: &mut XMbTables<'_>, ascending_cells: bool, bucket: usize, q_off_1: u32) {
    const K_COMPRESSION_FACTOR: usize = 2048;
    let old = t.hashtable[bucket];
    let empty = (old == 0) as u32;
    let idx = bucket >> t.pv_array_bts;
    if idx < t.pv_array.len() {
        t.pv_array[idx] |= empty << (bucket & PV_ARRAY_MASK);
    }
    let h = &mut t.helper_array[bucket / K_COMPRESSION_FACTOR];
    *h = h.saturating_add(1 - empty);
    t.next_pos[q_off_1 as usize] = old;
    t.hashtable[bucket] = q_off_1;
    if ascending_cells {
        t.words.push((bucket as u32, q_off_1));
    }
}

#[inline(always)]
fn apply<const H: u8>(t: &mut XMbTables<'_>, ascending_cells: bool, bucket: usize, q_off_1: u32) {
    if H == 4 {
        x_apply_bf(t, ascending_cells, bucket, q_off_1);
    } else {
        x_mb_apply(t, ascending_cells, bucket, q_off_1);
    }
}

/// The LOSAT_X_MBDELAY loop with the delay `D` and the prefetch hint `H` (D = 0: no queue).
#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn fill_delay<const D: usize, const H: u8>(
    queries_blastna: &[&[u8]],
    query_offsets: &[i32],
    query_masks: &[Vec<MaskedInterval>],
    word_length: usize,
    lut_word_length: usize,
    kmer_mask: u64,
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
    ascending_cells: bool,
    t: &mut XMbTables<'_>,
) {
    let mut x_ring = [(0u32, 0u32); D];
    let mut x_ring_len = 0usize;
    for_each_word(
        queries_blastna,
        query_offsets,
        query_masks,
        word_length,
        lut_word_length,
        kmer_mask,
        db_word_counts,
        max_db_word_count,
        t,
        |t, bucket, q_off_1| {
            if D == 0 {
                apply::<H>(t, ascending_cells, bucket, q_off_1);
                return;
            }
            prefetch::<H>(t.hashtable.as_ptr().wrapping_add(bucket));
            let slot = x_ring_len % D;
            if x_ring_len >= D {
                let (b, q) = x_ring[slot];
                apply::<H>(t, ascending_cells, b as usize, q);
            }
            x_ring[slot] = (bucket as u32, q_off_1);
            x_ring_len += 1;
        },
    );
    if D > 0 {
        let pending = x_ring_len.min(D);
        let first = x_ring_len - pending;
        for k in first..x_ring_len {
            let (b, q) = x_ring[k % D];
            apply::<H>(t, ascending_cells, b as usize, q);
        }
    }
}

/// Runs the probe's fill shape `(depth, hint)`; returns false for a shape that is not built.
#[allow(clippy::too_many_arguments)]
pub(super) fn fill(
    shape: (usize, u8),
    queries_blastna: &[&[u8]],
    query_offsets: &[i32],
    query_masks: &[Vec<MaskedInterval>],
    word_length: usize,
    lut_word_length: usize,
    kmer_mask: u64,
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
    ascending_cells: bool,
    t: &mut XMbTables<'_>,
) -> bool {
    macro_rules! go {
        ($d:literal, $h:literal) => {
            fill_delay::<$d, $h>(
                queries_blastna,
                query_offsets,
                query_masks,
                word_length,
                lut_word_length,
                kmer_mask,
                db_word_counts,
                max_db_word_count,
                ascending_cells,
                t,
            )
        };
    }
    match shape {
        (0, 4) => go!(0, 4),
        (0, _) => go!(0, 0),
        (1, 0) => go!(1, 0),
        (2, 0) => go!(2, 0),
        (4, 0) => go!(4, 0),
        (8, 0) => go!(8, 0),
        (16, 0) => go!(16, 0),
        (32, 0) => go!(32, 0),
        (64, 0) => go!(64, 0),
        (128, 0) => go!(128, 0),
        (16, 2) => go!(16, 2),
        (32, 2) => go!(32, 2),
        (64, 2) => go!(64, 2),
        (16, 3) => go!(16, 3),
        (32, 3) => go!(32, 3),
        (16, 4) => go!(16, 4),
        (32, 4) => go!(32, 4),
        (64, 4) => go!(64, 4),
        _ => return false,
    }
    true
}

/// `dump=<path>`: the fill's words, before the fill (no table is written).
#[allow(clippy::too_many_arguments)]
#[inline(never)]
pub(super) fn dump(
    queries_blastna: &[&[u8]],
    query_offsets: &[i32],
    query_masks: &[Vec<MaskedInterval>],
    word_length: usize,
    lut_word_length: usize,
    kmer_mask: u64,
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
    t: &mut XMbTables<'_>,
) {
    let Some(path) = probe().dump.as_ref() else {
        return;
    };
    let (tp, amb) = (*t.total_positions, *t.ambiguous_skipped);
    let mut out: Vec<u8> = Vec::new();
    for_each_word(
        queries_blastna,
        query_offsets,
        query_masks,
        word_length,
        lut_word_length,
        kmer_mask,
        db_word_counts,
        max_db_word_count,
        t,
        |_, bucket, q_off_1| {
            out.extend_from_slice(&(bucket as u32).to_le_bytes());
            out.extend_from_slice(&q_off_1.to_le_bytes());
        },
    );
    // the fill counts the words again
    *t.total_positions = tp;
    *t.ambiguous_skipped = amb;
    match std::fs::write(path, &out) {
        Ok(()) => eprintln!("[R10] dump words={} path={path}", out.len() / 8),
        Err(e) => eprintln!("[R10] dump failed {path}: {e}"),
    }
}
