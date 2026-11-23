//! dustmasker's masks of one sequence (`GetDustMasks_SkipNs`): runs of `N` are masked as they
//! are and the pieces between them go through one `CSymDustMasker` (`utils::dust`).

use crate::utils::dust::{DustMasker, NcbiLfgRandom};

/// A masked interval: the first and the last position (`CSymDustMasker::TMaskedInterval`,
/// two `TSeqPos`).
pub(super) type Interval = (u32, u32);

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:184-215
/// `s_FindSegmentWithLongNs` (with `max_ns` and `seq`): the runs of `N` that start the sequence,
/// end it, or are longer than `max_ns` (dustmasker's `-window` as given, as an unsigned 32-bit
/// number). NCBI counts each letter in `Uint4`s: a run starts at its index modulo 2^32 and is as
/// long as its length modulo 2^32, and a run whose length is a multiple of 2^32 counts as none.
///
/// The letters are read 16 at a time as a mask of their `N`s (`n_mask`). A block that leaves the
/// state as it was (no `N` after a letter that is not `N`, or only `N` inside a run) is passed
/// with that one test. In the others, the run that comes into the block is closed and the one
/// that goes on after it is opened, both at the edges of the mask. The runs that start and end
/// inside the block are visited only when one of them is kept (`inner_runs`): one longer than
/// `max_ns` (`RunTest`) or one at position 0.
fn long_n_runs(max_ns: u32, seq: &[u8]) -> Vec<Interval> {
    let mut runs = Vec::new();
    let close = |runs: &mut Vec<Interval>, start: usize, end: usize| {
        let (pos, ns) = (start as u32, (end - start) as u32);
        if ns > 0 && (ns > max_ns || pos == 0) {
            runs.push((pos, pos.wrapping_add(ns).wrapping_sub(1)));
        }
    };
    let longer = RunTest::new(max_ns);
    // Whether the letters so far end in a run of `N`, and the index where that run starts.
    let mut in_run = false;
    let mut open = 0;
    let (blocks, tail) = seq.as_chunks::<16>();
    for (index, block) in blocks.iter().enumerate() {
        let (low, high) = n_mask(block);
        if !in_run && (low | high) == 0 {
            continue;
        }
        let mask = gather(low) | gather(high) << 8;
        if in_run && mask == 0xffff {
            continue;
        }
        let base = index * 16;
        let others = !mask & 0xffff;
        let mut inner = mask;
        if in_run {
            let end = others.trailing_zeros();
            close(&mut runs, open, base + end as usize);
            inner &= !((1 << end) - 1);
        }
        in_run = mask & 0x8000 != 0;
        if in_run {
            let start = 32 - others.leading_zeros();
            open = base + start as usize;
            inner &= (1 << start) - 1;
        }
        let base = base as u32;
        if longer.holds(inner) != 0 || (base == 0 && inner & 1 != 0) {
            inner_runs(&mut runs, inner, base, max_ns);
        }
    }
    let base = seq.len() - tail.len();
    for (offset, &letter) in tail.iter().enumerate() {
        if letter == b'N' {
            if !in_run {
                (in_run, open) = (true, base + offset);
            }
        } else if in_run {
            close(&mut runs, open, base + offset);
            in_run = false;
        }
    }
    if in_run {
        let (pos, ns) = (open as u32, (seq.len() - open) as u32);
        if ns > 0 {
            runs.push((pos, pos.wrapping_add(ns).wrapping_sub(1)));
        }
    }
    runs
}

/// The high bit of each byte of `word` that is zero; no other bit is set.
fn zero_bytes(word: u64) -> u64 {
    const LOW7: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    !(((word & LOW7) + LOW7) | word | LOW7)
}

/// The `N`s of `block` as the high bits of its bytes, 8 letters in each word.
fn n_mask(block: &[u8; 16]) -> (u64, u64) {
    const NS: u64 = 0x4e4e_4e4e_4e4e_4e4e;
    let (low, high) = block.split_at(8);
    let word = |bytes: &[u8]| u64::from_le_bytes(bytes.try_into().expect("8 letters")) ^ NS;
    (zero_bytes(word(low)), zero_bytes(word(high)))
}

/// Bit k of the result is the high bit of byte k of `bits` (`bits` as `zero_bytes` gives).
fn gather(bits: u64) -> u32 {
    ((bits >> 7).wrapping_mul(0x0102_0408_1020_4080) >> 56) as u32
}

/// Whether a 16-bit mask holds a run of more than `max_ns` set bits: `holds(x)` has bit k set
/// when bits k to k + `max_ns` of `x` are all set (never when `max_ns` is 16 or more). Each shift
/// at most doubles the length of the runs it tests, so 5 shifts reach 17.
struct RunTest {
    shifts: [u32; 5],
}

impl RunTest {
    fn new(max_ns: u32) -> Self {
        let wanted = max_ns.saturating_add(1).min(17);
        let mut shifts = [0; 5];
        let mut covered = 1;
        for shift in &mut shifts {
            *shift = covered.min(wanted - covered);
            covered += *shift;
        }
        Self { shifts }
    }

    fn holds(&self, mut x: u32) -> u32 {
        for &shift in &self.shifts {
            x &= x >> shift;
        }
        x
    }
}

/// The kept runs of `inner`, a block's mask of `N`s without the runs that cross its edges, the
/// block starting at `base` (modulo 2^32): the k-th start of a run with the k-th end. A function
/// of its own, so that the loop over the blocks keeps its registers.
#[inline(never)]
fn inner_runs(runs: &mut Vec<Interval>, inner: u32, base: u32, max_ns: u32) {
    let mut starts = inner & !(inner << 1);
    let mut ends = (inner << 1) & !inner;
    while starts != 0 {
        let (start, end) = (starts.trailing_zeros(), ends.trailing_zeros());
        starts &= starts - 1;
        ends &= ends - 1;
        let pos = base.wrapping_add(start);
        if end - start > max_ns || pos == 0 {
            runs.push((pos, pos.wrapping_add(end - start - 1)));
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:212-225
/// `s_InsertMerge` (of `list`, `mask` and `linker`): `mask` extends the last interval when it
/// starts exactly `linker` (the `-linker` as given, unsigned 32-bit, the sum wrapping as NCBI's)
/// after its end.
fn insert_merge(list: &mut Vec<Interval>, mask: Interval, linker: u32) {
    if let Some(last) = list.last_mut() {
        if last.1.wrapping_add(linker) == mask.0 {
            last.1 = mask.1;
            return;
        }
    }
    list.push(mask);
}

/// `duster(seq, start, stop)` as intervals.
fn dust(
    duster: &DustMasker,
    seq: &[u8],
    start: u32,
    stop: u32,
    rng: &mut NcbiLfgRandom,
) -> Vec<Interval> {
    duster
        .mask_range(seq, start as usize, stop as usize, rng)
        .into_iter()
        .map(|interval| (interval.start as u32, (interval.end - 1) as u32))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:222-276
/// `GetDustMasks_SkipNs` (of `seq`, `level`, `window` and `linker`): the masks of the IUPAC letters
/// `seq` (not empty). The masker clamps `level`, `window` and `linker` to its ranges
/// (`DustMasker::new`, as `CSymDustMasker`'s constructor), while the runs of `N` and the merges
/// take `window` and `linker` as given (`Uint4` of the integer arguments, so `-1` is
/// 4,294,967,295). One masker, and so one sequence of random numbers for `N`, serves every piece of
/// the sequence.
pub(super) fn dust_masks_skip_ns(
    seq: &[u8],
    level: u32,
    window: u32,
    linker: u32,
) -> Vec<Interval> {
    let duster = DustMasker::new(level, window as usize, linker as usize);
    let mut rng = NcbiLfgRandom::new();
    let last = (seq.len() - 1) as u32;
    let n_runs = long_n_runs(window, seq);
    if n_runs.is_empty() {
        return dust(&duster, seq, 0, last, &mut rng);
    }
    let mut rv: Vec<Interval> = Vec::new();
    let mut seq_start: u32 = 0;
    for &run in &n_runs {
        if run.0 == 0 {
            seq_start = run.1.wrapping_add(1);
            rv.push(run);
            continue;
        }
        let s_mask = dust(&duster, seq, seq_start, run.0 - 1, &mut rng);
        if let Some((&first, rest)) = s_mask.split_first() {
            insert_merge(&mut rv, first, linker);
            rv.extend_from_slice(rest);
            insert_merge(&mut rv, run, linker);
        } else {
            rv.push(run);
        }
        seq_start = run.1.wrapping_add(1);
    }
    if (seq_start as usize) < seq.len() {
        let s_mask = dust(&duster, seq, seq_start, last, &mut rng);
        if let Some((&first, rest)) = s_mask.split_first() {
            insert_merge(&mut rv, first, linker);
            rv.extend_from_slice(rest);
        }
    }
    rv
}

#[cfg(test)]
#[path = "masker_tests.rs"]
mod tests;
