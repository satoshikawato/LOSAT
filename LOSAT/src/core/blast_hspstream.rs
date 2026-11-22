//! The close of the preliminary stage's HSP stream with composition-based statistics and
//! `ADAPTIVE_CBS` (NCBI core/blast_hspstream.c `BlastHSPCBSStreamClose`): each query's
//! preliminary hit list sorted by E-value and cut where the E-values pass a limit, before the
//! traceback. BLASTP, PSIBLAST (whose traceback never takes it), RPSBLAST, BLASTX and TBLASTN
//! call it, each through a `cbs_close_trim` of its hit list type that sorts the list its own way
//! between `cbs_trim_sorts` and `cbs_trim_count`. The environment variable is read here only
//! (`adaptive_cbs_environment`), also for the preliminary hit list size
//! (`blastn::hsp::get_prelim_hitlist_size`).

use crate::utils::glibc_fma_libm::GlibcLibm;

/// Whether the environment variable `ADAPTIVE_CBS` is set (to any value, the empty one too):
/// with composition-based statistics NCBI then keeps a larger preliminary hit list
/// (`blastn::hsp::get_prelim_hitlist_size`) and closes the stream with
/// `BlastHSPCBSStreamClose` (`adaptive_cbs_close`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:47
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
pub(crate) fn adaptive_cbs_close(composition_based_stats: bool, hitlist_size: usize) -> bool {
    composition_based_stats && hitlist_size < 1000 && adaptive_cbs_environment()
}

/// The test of `BlastHSPCBSStreamClose` before it sorts one query's preliminary hit list: the
/// list holds more than `REF_HIT_NUM` plus 700 HSP lists (`hitlist_size` is the traceback's hit
/// list size, below 1000; `cbs_trim_count` follows).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:827-838
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
/// `mag` takes a `log10` (a `double`) converted to `int` toward zero, and `mag` times `pct` over
/// 100 and `mag` over 2 are C's `int` divisions (toward zero), as Rust's `as i32` and `/` are.
/// `max_index` is taken before the sort, whose purge could only shorten a list that held empty HSP
/// lists (the preliminary stage keeps none); the loop also stops at the sorted count, where NCBI
/// would read past it.
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
        mag = ref_evalue.glibc_log10() as i32;
    }
    if mag < -20 {
        mag = (mag * 90 / 100).max(mag + 10);
    } else {
        mag /= 2;
    }
    let evalue_limit = if mag >= 0 {
        ref_evalue * 3.0
    } else {
        9.9 * 10f64.glibc_pow(f64::from(mag))
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
#[path = "blast_hspstream_tests.rs"]
mod tests;
