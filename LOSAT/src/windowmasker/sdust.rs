//! `-dust T`: `CSDustMasker` (the symmetric DUST masker of windowmasker, run over the parts of
//! a sequence that WindowMasker did not mask) and `CSeqMasker::MergeMaskInfo`.

use crate::utils::dust::{DustMasker, NcbiLfgRandom};

/// `CSDustMasker`: one `CSymDustMasker` (window 64, the level of `-dust_level`, linker 1) for
/// the whole run, so the random numbers that resolve `N` run on from one sequence to the next.
pub(super) struct SDustMasker {
    duster: DustMasker,
    rng: NcbiLfgRandom,
}

/// `CWinMaskConfig`'s `dust_window` (64) and `dust_linker` (1).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:358-361
const DUST_WINDOW: u32 = 64;
const DUST_LINKER: usize = 1;

impl SDustMasker {
    /// `new CSDustMasker( aConfig.DustWindow(), aConfig.DustLevel(), aConfig.DustLinker() )`:
    /// the level is the `Uint4` of `-dust_level` (a negative value wraps, and `CSymDustMasker`
    /// takes 20 for a level outside 2..64, `DustMasker::new`).
    ///
    /// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_sdust_masker.cpp:67-69
    pub(super) fn new(level: i32) -> Self {
        SDustMasker {
            duster: DustMasker::new(level as u32, DUST_WINDOW as usize, DUST_LINKER),
            rng: NcbiLfgRandom::new(),
        }
    }

    /// `(*duster)( data, *mask_info )`: the DUST intervals of `data` (IUPACna) outside the
    /// WindowMasker intervals `exclude` longer than the window, both ends included.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_sdust_masker.cpp:76-116
    /// `CSeqVector_CI::operator+` sets the position through `SetPos`, which stops at the
    /// sequence's end: `end_it` is at most `data.end()`. The duster's `stop` is the last
    /// position it reads (`CSymDustMasker` clamps it to the last letter).
    pub(super) fn mask(&mut self, data: &[u8], exclude: &[(u32, u32)]) -> Vec<(u32, u32)> {
        let len = data.len() as u64;
        let window = u64::from(DUST_WINDOW);
        let mut result = Vec::new();
        let mut e = exclude.iter().peekable();
        let mut current = 0u64;
        if let Some(&&(first, second)) = e.peek() {
            if first == 0 && u64::from(second) + 1 > window {
                current = u64::from(second) - window + 2;
                e.next();
            }
        }
        loop {
            let start = current;
            while e
                .peek()
                .is_some_and(|&&(first, second)| u64::from(second - first) < window)
            {
                e.next();
            }
            let end = match e.peek() {
                None => len,
                Some(&&(first, _)) => (u64::from(first) + window).min(len),
            };
            result.extend(
                self.duster
                    .mask_range(data, start as usize, end as usize, &mut self.rng)
                    .into_iter()
                    .map(|interval| (interval.start as u32, (interval.end - 1) as u32)),
            );
            if let Some(&(_, second)) = e.next() {
                current = u64::from(second) - window + 2;
            }
            if end == len {
                break;
            }
        }
        result
    }
}

/// `CSeqMasker::MergeMaskInfo( dest, src )`: the intervals of `dest` and `src` (each in order),
/// those that overlap or touch merged, into `dest`.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:507-554
/// (`seg.second + 1` is a `TSeqPos` sum; positions are below 2^32 - 1, so it does not wrap.)
pub(super) fn merge_mask_info(dest: &mut Vec<(u32, u32)>, src: &[(u32, u32)]) {
    if src.is_empty() {
        return;
    }
    let (mut si, mut di) = (0, 0);
    let mut seg = if dest.first().is_some_and(|d| d.0 < src[0].0) {
        di += 1;
        dest[0]
    } else {
        si += 1;
        src[0]
    };
    let mut res = Vec::with_capacity(dest.len() + src.len());
    loop {
        let next = match (src.get(si), dest.get(di)) {
            (Some(&s), Some(&d)) => {
                if s.0 < d.0 {
                    si += 1;
                    s
                } else {
                    di += 1;
                    d
                }
            }
            (Some(&s), None) => {
                si += 1;
                s
            }
            (None, Some(&d)) => {
                di += 1;
                d
            }
            (None, None) => break,
        };
        if u64::from(seg.1) + 1 < u64::from(next.0) {
            res.push(seg);
            seg = next;
        } else if seg.1 < next.1 {
            seg.1 = next.1;
        }
    }
    res.push(seg);
    *dest = res;
}

#[cfg(test)]
#[path = "sdust_tests.rs"]
mod tests;
