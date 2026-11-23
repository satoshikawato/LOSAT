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
/// The search goes 16 letters at a time (`first_from`).
fn long_n_runs(max_ns: u32, seq: &[u8]) -> Vec<Interval> {
    let mut runs = Vec::new();
    let mut from = 0;
    while let Some(start) = first_from(seq, from, |letter| letter == b'N') {
        let end = first_from(seq, start, |letter| letter != b'N').unwrap_or(seq.len());
        let (pos, ns) = (start as u32, (end - start) as u32);
        if ns > 0 && (end == seq.len() || ns > max_ns || pos == 0) {
            runs.push((pos, pos.wrapping_add(ns).wrapping_sub(1)));
        }
        from = end;
    }
    runs
}

/// The index of the first letter of `seq` from `from` on for which `wanted` holds. Each block of
/// 16 letters is tested whole first (a test without a branch per letter).
fn first_from(seq: &[u8], from: usize, wanted: impl Fn(u8) -> bool) -> Option<usize> {
    let rest = &seq[from..];
    let (blocks, _) = rest.as_chunks::<16>();
    let at = blocks
        .iter()
        .position(|block| {
            block
                .iter()
                .fold(false, |any, &letter| any | wanted(letter))
        })
        .map_or(blocks.len() * 16, |block| block * 16);
    rest[at..]
        .iter()
        .position(|&letter| wanted(letter))
        .map(|off| from + at + off)
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
