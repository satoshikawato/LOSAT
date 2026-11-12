//! dustmasker's masks of one sequence (`GetDustMasks_SkipNs`): runs of `N` are masked as they
//! are and the pieces between them go through one `CSymDustMasker` (`utils::dust`).

use crate::utils::dust::{DustMasker, NcbiLfgRandom};

/// A masked interval: the first and the last position (`CSymDustMasker::TMaskedInterval`,
/// two `TSeqPos`).
pub(super) type Interval = (u32, u32);

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:184-210
// ```c++
// CSymDustMasker::TMaskList s_FindSegmentWithLongNs(const unsigned int MAX_Ns, objects::CSeqVector & seq)
// {
// 	//Always trim Ns at start and end of seq
// 	CSymDustMasker::TMaskList  NsRange;
//     unsigned int pos = 0;
// 	unsigned int Ns = 0;
//     for(objects::CSeqVector_CI itr=seq.begin();itr!=seq.end(); ++itr) {
//         if ((*itr) == 78) {
//             Ns++;
//         }
//         else {
//             if (Ns > 0) {
//             	if((Ns > MAX_Ns ) || (pos == 0)) {
//             		CSymDustMasker::TMaskedInterval r(pos, pos+Ns-1);
//             		NsRange.push_back(r);
//             	}
//                 pos += Ns;
//                 Ns = 0;
//             }
//            	pos++;
//         }
//     }
//     if (Ns > 0) {
//     	CSymDustMasker::TMaskedInterval r(pos, pos+Ns-1);
//     	NsRange.push_back(r);
//     }
//     return NsRange;
// }
// ```
/// `s_FindSegmentWithLongNs(max_ns, seq)`: the runs of `N` that start the sequence, end it,
/// or are longer than `max_ns` (dustmasker's `-window` as given, as an unsigned 32-bit number).
fn long_n_runs(max_ns: u32, seq: &[u8]) -> Vec<Interval> {
    let mut runs = Vec::new();
    let mut pos: u32 = 0;
    let mut ns: u32 = 0;
    for &letter in seq {
        if letter == b'N' {
            ns = ns.wrapping_add(1);
        } else {
            if ns > 0 {
                if ns > max_ns || pos == 0 {
                    runs.push((pos, pos.wrapping_add(ns).wrapping_sub(1)));
                }
                pos = pos.wrapping_add(ns);
                ns = 0;
            }
            pos = pos.wrapping_add(1);
        }
    }
    if ns > 0 {
        runs.push((pos, pos.wrapping_add(ns).wrapping_sub(1)));
    }
    runs
}

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:212-220
// ```c++
// void s_InsertMerge(CSymDustMasker::TMaskList & list, CSymDustMasker::TMaskedInterval & new_mask, Uint4 linker)
// {
// 	if ((!list.empty()) && (list.back().second + linker == new_mask.first)) {
// 		list.back().second = new_mask.second;
// 		return;
// 	}
// 	list.push_back(new_mask);
// }
// ```
/// `s_InsertMerge(list, mask, linker)`: `mask` extends the last interval when it starts exactly
/// `linker` (the `-linker` as given, unsigned 32-bit, the sum wrapping as NCBI's) after its end.
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
// ```c++
// std::unique_ptr< CSymDustMasker::TMaskList >
// GetDustMasks_SkipNs(objects::CSeqVector & seq, Uint4 level, Uint4 window, Uint4 linker)
// {
//     CSymDustMasker duster(level, window, linker);
//     CSymDustMasker::TMaskList NsRange = s_FindSegmentWithLongNs(window, seq);
//     if(NsRange.empty()){
//         return duster(seq);
//     }
//
//     std::unique_ptr< CSymDustMasker::TMaskList > rv(new CSymDustMasker::TMaskList);
//     TSeqPos seq_start =0;
//     NON_CONST_ITERATE(CSymDustMasker::TMaskList, itr, NsRange) {
//     	if(itr->first == 0) {
//     		seq_start = itr->second + 1;
//     		rv->push_back(*itr);
//     		continue;
//     	}
//     	else {
//     		std::unique_ptr< CSymDustMasker::TMaskList > s_mask = duster(seq, seq_start, itr->first -1);
//     		if(s_mask->size() > 0) {
//     			s_InsertMerge(*rv, s_mask->front(), linker);
//     			if( s_mask->size() > 1) {
//     				rv->insert(rv->end(), ++(s_mask->begin()), s_mask->end());
//     			}
//     			s_InsertMerge(*rv, *itr, linker);
//     		}
//     		else {
//     			rv->push_back(*itr);
//     		}
//    			seq_start = itr->second + 1;
//     	}
//     }
//
//     if(seq_start < seq.size()){
//    		std::unique_ptr< CSymDustMasker::TMaskList > s_mask = duster(seq, seq_start, seq.size() -1);
//    		if(s_mask->size() > 0) {
//    			s_InsertMerge(*rv, s_mask->front(), linker);
//    			if( s_mask->size() > 1) {
//    				rv->insert(rv->end(), ++(s_mask->begin()), s_mask->end());
//    			}
//    		}
//     }
//
//     return rv;
// }
// ```
/// `GetDustMasks_SkipNs(seq, level, window, linker)`: the masks of the IUPAC letters `seq`
/// (not empty). The masker clamps `level`, `window` and `linker` to its ranges
/// (`DustMasker::new`, as `CSymDustMasker`'s constructor), while the runs of `N` and the merges
/// take `window` and `linker` as given (`Uint4` of the integer arguments, so `-1` is
/// 4,294,967,295). One masker, and so one sequence of random numbers for `N`, serves every
/// piece of the sequence.
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
mod tests {
    use super::*;

    #[test]
    fn runs_of_n_at_the_ends_always_count_and_inner_ones_only_when_longer_than_max() {
        assert_eq!(
            long_n_runs(2, b"NNACNNTNNNG"),
            vec![(0, 1), (7, 9)],
            "the run of 2 inside is not longer than 2"
        );
        assert_eq!(long_n_runs(64, b"ACGTN"), vec![(4, 4)]);
        assert_eq!(long_n_runs(u32::MAX, b"NNNN"), vec![(0, 3)]);
        assert!(long_n_runs(0, b"ACGT").is_empty());
    }

    #[test]
    fn merges_take_the_linker_as_an_unsigned_sum() {
        let mut list = vec![(0, 9)];
        insert_merge(&mut list, (10, 20), 1);
        assert_eq!(list, vec![(0, 20)]);
        // -1 as Uint4: 20 + 4294967295 wraps to 19.
        insert_merge(&mut list, (19, 30), u32::MAX);
        assert_eq!(list, vec![(0, 30)]);
        insert_merge(&mut list, (40, 50), 1);
        assert_eq!(list, vec![(0, 30), (40, 50)]);
    }
}
