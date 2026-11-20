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
/// ```c++
///       dust_window( 64 ),
///       dust_level( app_type == eGenerateMasksWithDuster ? args["dust_level"].AsInteger() : 0 ),
///       dust_linker( 1 ),
/// ```
const DUST_WINDOW: u32 = 64;
const DUST_LINKER: usize = 1;

impl SDustMasker {
    /// `new CSDustMasker( aConfig.DustWindow(), aConfig.DustLevel(), aConfig.DustLinker() )`:
    /// the level is the `Uint4` of `-dust_level` (a negative value wraps, and `CSymDustMasker`
    /// takes 20 for a level outside 2..64, `DustMasker::new`).
    ///
    /// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_sdust_masker.cpp:67-69
    /// ```c++
    /// CSDustMasker::CSDustMasker( Uint4 arg_window, Uint4 arg_level, Uint4 arg_linker )
    ///   : window( arg_window ), duster_( arg_level, arg_window, arg_linker )
    /// {}
    /// ```
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
    /// ```c++
    ///     TMaskList * res( new TMaskList );
    ///     TMaskList::const_iterator e_it = exclude_ranges.begin();
    ///     TMaskList::const_iterator e_end = exclude_ranges.end();
    ///     CSeqVector::const_iterator start_it = data.begin();
    ///     CSeqVector::const_iterator current_it = data.begin();
    ///     CSeqVector::const_iterator end_it = data.end();
    ///
    ///     if( e_it != e_end && e_it->first == 0 && e_it->second + 1 > window )
    ///     {
    ///         current_it = start_it + e_it->second - window + 2;
    ///         ++e_it;
    ///     }
    ///
    ///     do
    ///     {
    ///         Uint4 start_offset = current_it - start_it;
    ///
    ///         while( e_it != e_end && e_it->second - e_it->first + 1 <= window )
    ///             ++e_it;
    ///
    ///         end_it = (e_it == e_end) ? data.end()
    ///                                    : start_it + e_it->first + window;
    ///         Uint4 stop_offset = end_it - start_it;
    ///
    ///         // Now dust.
    ///         std::unique_ptr< TMaskList > result = duster_(
    ///             data, start_offset, stop_offset );
    ///
    ///         res->insert( res->end(), result->begin(), result->end() );
    ///
    ///         if( e_it != e_end )
    ///         {
    ///             current_it = start_it + e_it->second - window + 2;
    ///             ++e_it;
    ///         }
    ///     } while( end_it != data.end() );
    /// ```
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
                .is_some_and(|&&(first, second)| u64::from(second - first) + 1 <= window)
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
/// ```c++
///     if( src->empty() )
///         return;
///
///     TMaskList::const_iterator si( src->begin() );
///     TMaskList::const_iterator send( src->end() );
///     TMaskList::iterator di( dest->begin() );
///     TMaskList::iterator dend( dest->end() );
///     TMaskList res;
///     TMaskedInterval seg;
///     TMaskedInterval next_seg;
///
///     if( di != dend && di->first < si->first )
///         seg = *(di++);
///     else seg = *(si++);
///
///     while( true )
///     {
///         if( si != send ) {
///             if( di != dend ) {
///                 if( si->first < di->first ) {
///                     next_seg = *(si++);
///                 } else {
///                     next_seg = *(di++);
///                 }
///             } else {
///                 next_seg = *(si++);
///             }
///         } else if( di != dend ) {
///             next_seg = *(di++);
///         } else {
///             break;
///         }
///
///         if( seg.second + 1 < next_seg.first ) {
///             res.push_back( seg );
///             seg = next_seg;
///         }
///         else if( seg.second < next_seg.second ) {
///             seg.second = next_seg.second;
///         }
///     }
///
///     res.push_back( seg );
///     dest->swap(res);
/// ```
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
mod tests {
    use super::*;

    #[test]
    fn merging_joins_overlapping_and_touching_intervals() {
        let mut dest = vec![(0, 9), (30, 39)];
        merge_mask_info(&mut dest, &[(10, 12), (20, 25), (35, 50)]);
        assert_eq!(dest, vec![(0, 12), (20, 25), (30, 50)]);
        let mut dest = vec![];
        merge_mask_info(&mut dest, &[(1, 2)]);
        assert_eq!(dest, vec![(1, 2)]);
        let mut dest = vec![(1, 2)];
        merge_mask_info(&mut dest, &[]);
        assert_eq!(dest, vec![(1, 2)]);
    }

    #[test]
    fn a_sequence_without_masks_is_dusted_whole() {
        let data = b"ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC";
        let mut masker = SDustMasker::new(20);
        let whole = masker.mask(data, &[]);
        let direct = DustMasker::new(20, 64, 1)
            .mask_range(data, 0, data.len(), &mut NcbiLfgRandom::new())
            .into_iter()
            .map(|interval| (interval.start as u32, (interval.end - 1) as u32))
            .collect::<Vec<_>>();
        assert_eq!(whole, direct);
        assert!(!whole.is_empty());
    }
}
