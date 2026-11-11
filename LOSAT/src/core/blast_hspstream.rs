//! The close of the preliminary stage's HSP stream with composition-based statistics and
//! `ADAPTIVE_CBS` (NCBI core/blast_hspstream.c `BlastHSPCBSStreamClose`): each query's
//! preliminary hit list sorted by E-value and cut where the E-values pass a limit, before the
//! traceback. BLASTP, PSIBLAST (whose traceback never takes it), RPSBLAST, BLASTX and TBLASTN
//! call it, each through a `cbs_close_trim` of its hit list type that sorts the list its own way
//! between `cbs_trim_sorts` and `cbs_trim_count`. The environment variable is read here only
//! (`adaptive_cbs_environment`), also for the preliminary hit list size
//! (`blastn::hsp::get_prelim_hitlist_size`).

/// Whether the environment variable `ADAPTIVE_CBS` is set (to any value, the empty one too):
/// with composition-based statistics NCBI then keeps a larger preliminary hit list
/// (`blastn::hsp::get_prelim_hitlist_size`) and closes the stream with
/// `BlastHSPCBSStreamClose` (`adaptive_cbs_close`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:51
/// ```c
///     char * ADAPTIVE_CBS_ENV = getenv("ADAPTIVE_CBS");
/// ```
pub(crate) fn adaptive_cbs_environment() -> bool {
    std::env::var_os("ADAPTIVE_CBS").is_some()
}

/// Whether the traceback starts with `BlastHSPCBSStreamClose` (each query's hit list's
/// `cbs_close_trim`) in place of `BlastHSPStreamClose`: `ADAPTIVE_CBS` is set, composition-based
/// statistics are on and the traceback's hit list size is below 1000. PSIBLAST's traceback
/// takes the preliminary size (traceback_stage.cpp:216-223), 1500 or more here, so it never
/// trims.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1825,1839-1846
/// ```c
///     char * ADAPTIVE_CBS_ENV = getenv("ADAPTIVE_CBS");
/// ...
///     /* Prohibit any subsequent writing to the HSP stream. */
///     if((ADAPTIVE_CBS_ENV != NULL) && (ext_options->compositionBasedStats) && (hit_options->hitlist_size < 1000)) {
///     	BlastHSPCBSStreamClose(hsp_stream, hit_options->hitlist_size);
///
///     }
///     else {
///     	BlastHSPStreamClose(hsp_stream);
///     }
/// ```
pub(crate) fn adaptive_cbs_close(composition_based_stats: bool, hitlist_size: usize) -> bool {
    composition_based_stats && hitlist_size < 1000 && adaptive_cbs_environment()
}

/// The test of `BlastHSPCBSStreamClose` before it sorts one query's preliminary hit list: the
/// list holds more than `REF_HIT_NUM + 700` HSP lists (`hitlist_size` is the traceback's hit
/// list size, below 1000; `cbs_trim_count` follows).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:827-835
/// ```c
///       hit_list = results->hitlist_array[index];
///       if (hit_list) {
///     	  const int REF_HIT_NUM = MAX(500, hitlist_size);
///     	  const int MIN_BUF_SZ = REF_HIT_NUM + 600;
///     	  if (MIN_BUF_SZ + 100 < hit_list->hsplist_count) {
///         	  int max_index = hit_list->hsplist_count -1;
///         	  double best_evalue = 0, evalue_limit = 0;
///         	  int mag = -180, pct = 90, i = 0;
///         	  Blast_HitListSortByEvalue(hit_list);
/// ```
pub(crate) fn cbs_trim_sorts(hsplist_count: usize, hitlist_size: usize) -> bool {
    hitlist_size.max(500) + 700 < hsplist_count
}

/// The rest of `BlastHSPCBSStreamClose` for one query's preliminary hit list, after
/// `Blast_HitListSortByEvalue`: the count the list is cut to (`s_TrimHitList`), or `None`.
/// `count` is the list's count before the sort (`cbs_trim_sorts` held for it), `sorted` its count
/// after the sort's purge of empty HSP lists, and `best_evalue(i)` the best E-value of its `i`-th
/// HSP list in the sorted order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:806-864
/// ```c
/// void s_TrimHitList(BlastHitList* hitlist, int count)
/// {
/// 	int old_count = hitlist->hsplist_count;
/// 	int index =0;
/// 	for (index = count; index < old_count; index++) {
/// 		hitlist->hsplist_array[index] = Blast_HSPListFree(hitlist->hsplist_array[index]);
/// 	}
/// 	hitlist->hsplist_count = count;
/// }
///
/// void BlastHSPCBSStreamClose(BlastHSPStream* hsp_stream, int hitlist_size)
/// {
/// ...
///         	  Blast_HitListSortByEvalue(hit_list);
///         	  best_evalue = hit_list->hsplist_array[REF_HIT_NUM]->best_evalue;
///         	  if (best_evalue != 0 ){
///         		  mag = log10(hit_list->hsplist_array[REF_HIT_NUM]->best_evalue);
///         	  }
///         	  if (mag < -20){
///         		  mag = MAX(mag*pct/100, mag +10);
///         	  }
///         	  else {
///         		 mag = mag/2;
///         	  }
///         	  evalue_limit = (mag >= 0)? best_evalue*3: 9.9* pow(10, mag);
///
///         	  i = MIN_BUF_SZ;
///         	  for(; i < max_index; i +=100) {
///         		  if((hit_list->hsplist_array[i]->best_evalue != 0) &&
///         		     (evalue_limit < hit_list->hsplist_array[i]->best_evalue)){
///         			  s_TrimHitList(hit_list, i);
///         			  break;
///         		  }
///         	  }
///           }
///       }
///    }
///
///    BlastHSPStreamClose(hsp_stream);
///    return;
///
/// }
/// ```
/// `mag = log10(...)` converts the `double` to `int` toward zero, and `mag*pct/100` and `mag/2`
/// are C's `int` divisions (toward zero), as Rust's `as i32` and `/` are. `max_index` is taken
/// before the sort, whose purge could only shorten a list that held empty HSP lists (the
/// preliminary stage keeps none); the loop also stops at the sorted count, where NCBI would read
/// past it.
pub(crate) fn cbs_trim_count(
    count: usize,
    sorted: usize,
    hitlist_size: usize,
    best_evalue: impl Fn(usize) -> f64,
) -> Option<usize> {
    let ref_hit_num = hitlist_size.max(500);
    let min_buf_sz = ref_hit_num + 600;
    if ref_hit_num >= sorted {
        return None;
    }
    let max_index = count - 1;
    let ref_evalue = best_evalue(ref_hit_num);
    let mut mag: i32 = -180;
    if ref_evalue != 0.0 {
        mag = ref_evalue.log10() as i32;
    }
    if mag < -20 {
        mag = (mag * 90 / 100).max(mag + 10);
    } else {
        mag /= 2;
    }
    let evalue_limit = if mag >= 0 {
        ref_evalue * 3.0
    } else {
        9.9 * 10f64.powf(f64::from(mag))
    };
    let mut i = min_buf_sz;
    while i < max_index && i < sorted {
        let evalue = best_evalue(i);
        if evalue != 0.0 && evalue_limit < evalue {
            return Some(i);
        }
        i += 100;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::{cbs_trim_count, cbs_trim_sorts};

    /// The count a hit list of the best E-values `evalues` (in E-value order) keeps after the
    /// trim, as a caller runs it: the test, the sort (here the order given) and the cut.
    fn trimmed(evalues: &[f64], hitlist_size: usize) -> usize {
        if !cbs_trim_sorts(evalues.len(), hitlist_size) {
            return evalues.len();
        }
        cbs_trim_count(evalues.len(), evalues.len(), hitlist_size, |i| evalues[i])
            .unwrap_or(evalues.len())
    }

    /// `count` E-values: `low` for the first `split`, then `high`.
    fn list(count: usize, split: usize, low: f64, high: f64) -> Vec<f64> {
        (0..count)
            .map(|i| if i < split { low } else { high })
            .collect()
    }

    // No NCBI unit test covers BlastHSPCBSStreamClose; the cases follow its arithmetic
    // (blast_hspstream.c:829-856) at the edges of each branch.
    #[test]
    fn the_list_is_cut_at_the_first_hundred_past_the_limit() {
        // 1,200 lists are not more than MIN_BUF_SZ + 100 = 1,200: nothing happens.
        assert!(!cbs_trim_sorts(1200, 500));
        assert_eq!(trimmed(&list(1200, 700, 1e-30, 1e-5), 500), 1200);
        // 1,201: E[500] = 1e-30, mag -30 -> max(-27, -20) = -20, limit 9.9e-20; E[1100] =
        // 1e-5 passes it: the list keeps 1,100.
        assert!(cbs_trim_sorts(1201, 500));
        assert_eq!(trimmed(&list(1201, 700, 1e-30, 1e-5), 500), 1100);
        // The first index past the limit is 1,250: the cut is at 1,300 (steps of 100), and
        // 1,300 is max_index of a 1,301-list, which the loop does not reach.
        assert_eq!(trimmed(&list(1400, 1250, 1e-30, 1e-5), 500), 1300);
        assert_eq!(trimmed(&list(1301, 1250, 1e-30, 1e-5), 500), 1301);
        // A larger hit list size moves both indices: REF 700, MIN_BUF_SZ 1,300.
        assert!(!cbs_trim_sorts(1400, 700));
        assert_eq!(trimmed(&list(1500, 900, 1e-30, 1e-5), 700), 1300);
        // A hit list size below 500 counts as 500.
        assert_eq!(trimmed(&list(1201, 700, 1e-30, 1e-5), 10), 1100);
    }

    #[test]
    fn the_limit_follows_ncbi_s_integer_magnitudes() {
        // E[500] = 1e-10: mag -10 -> -10/2 = -5, limit 9.9e-5; 1e-4 passes it, 9e-5 does not.
        assert_eq!(trimmed(&list(1300, 1000, 1e-10, 1e-4), 500), 1100);
        assert_eq!(trimmed(&list(1300, 1000, 1e-10, 9e-5), 500), 1300);
        // E[500] = 3e-25: log10 = -24.5 -> -24 (toward zero), -24*90/100 = -21 (toward zero),
        // max(-21, -14) = -14, limit 9.9e-14.
        assert_eq!(trimmed(&list(1300, 1000, 3e-25, 1e-13), 500), 1100);
        assert_eq!(trimmed(&list(1300, 1000, 3e-25, 9e-14), 500), 1300);
        // E[500] = 5: mag 0 -> 0, limit 3 * 5 = 15.
        assert_eq!(trimmed(&list(1300, 1000, 5.0, 16.0), 500), 1100);
        // E[500] = 0: mag stays -180 -> max(-162, -170) = -162; an E-value of 0 never cuts.
        assert_eq!(trimmed(&list(1300, 1300, 0.0, 0.0), 500), 1300);
        assert_eq!(trimmed(&list(1300, 1000, 0.0, 1e-160), 500), 1100);
    }

    #[test]
    fn a_list_shortened_by_the_purge_is_read_to_its_sorted_end() {
        // 1,300 lists of which 150 were empty: the loop reads indices 1,100 and below 1,150.
        let evalues = list(1150, 1000, 1e-30, 1e-5);
        assert_eq!(cbs_trim_count(1300, 1150, 500, |i| evalues[i]), Some(1100));
        let late = list(1150, 1120, 1e-30, 1e-5);
        assert_eq!(cbs_trim_count(1300, 1150, 500, |i| late[i]), None);
    }
}
