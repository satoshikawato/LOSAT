//! EXPERIMENT (LOSAT_X_BNSHRINKFREE; in the strict set: on unless it is `0`): the largest buffers a
//! query batch frees, shrunk to one element before they are freed: the CSR of a lookup table
//! (`x_dc_csr::Csr`: `offsets`, one element per cell, and `qs`, one per word; LOSAT_X_DCCSR's
//! discontiguous tables and LOSAT_X_BNCSR's contiguous ones) and the batch's BLASTNA query buffer
//! (`search_query_batch`'s `encoded_query_buffer`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1262,1376
//! ```c
//!    mb_lt->hashtable = (Int4*)calloc(mb_lt->hashsize, sizeof(Int4));
//! ...
//!    sfree(mb_lt->hashtable);
//! ```
//! NCBI allocates a batch's tables and query buffer for the batch and frees them with it, as LOSAT
//! does; how the memory goes back changes here, nothing else. No NCBI counterpart: memory only.
//!
//! What it changes. glibc maps a block of at least its mmap threshold (128 KiB at the start) and
//! unmaps it when it is freed, but a `free` of a mapped block also raises the threshold to the
//! block's size (up to 32 MiB), and the heap then serves every later allocation below it.
//! dc-megablast of 31,168 cDNAs at one thread freed a CSR of 16 + 36 MiB and a 10 MiB query buffer
//! after each of its 18 batches: after the first batch the threshold was 16 MiB, the heap took the
//! batches' other buffers, and its arena grew from 44.7 to 58.6 MiB with 13-20 MiB free in it,
//! which stayed resident (peak `RssAnon` 90.8 MiB against 80.8 for the memory floor of plan v4
//! line MF; with the threshold fixed by `GLIBC_TUNABLES` the same search peaked at about 79). With
//! the switch these buffers are shrunk to one element first (`free`: `shrink_to(1)`, which glibc's
//! `realloc` does by `mremap` for a mapped block, down to one page, without the `free` that moves
//! the threshold), then freed (the one-page block is below any threshold). The threshold stays
//! where the rest of the search puts it, and these buffers keep coming from fresh mappings, as at
//! the first batch. Keeping them for the next batch instead (the same address space, no page
//! faults) raised the peak by 10 MiB: LOSAT_X_BNLUTDROP frees a batch's tables before its
//! traceback, and a kept table stays resident through it.
//!
//! Why the outputs are the same. Only how the freed memory goes back changes; no element is read
//! after the shrink. No floating-point operation is involved.
//!
//! The shrink costs one `mremap` per buffer per batch. On a heap block (below the threshold, or
//! another allocator) `realloc` shortens it in place or copies one element, and the `free` follows
//! as before.
//!
//! Timing: with LOSAT_TIMING, `[TIMING] x_bn_shrink_free: buffers=` gives the number of buffers
//! shrunk before their free.
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: LOSAT_X_BNSHRINKFREE (read once).
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNSHRINKFREE"))
}

// No NCBI counterpart: the buffers shrunk before their free.
static SHRUNK: AtomicU64 = AtomicU64::new(0);

/// Frees `vector`, shrunk to one element first with the switch (module comment).
pub(crate) fn free<T>(mut vector: Vec<T>) {
    if on() && vector.capacity() > 1 {
        vector.clear();
        vector.shrink_to(1);
        SHRUNK.fetch_add(1, Ordering::Relaxed);
    }
    drop(vector);
}

/// PROBE (not for a pull request): `free` where bit `bit` of LOSAT_PROBE_SHRINKMASK is set
/// (unset: every bit), else a plain drop. 1 BLASTNA buffer, 2 CSR offsets (not kept), 4 CSR qs,
/// 8 chains' heads, 16 next_pos, 32 pv_array, 64 MB24 heads, 128 MB24 chain, 256 second template.
pub(crate) fn free_as<T>(bit: u32, vector: Vec<T>) {
    static MASK: OnceLock<u32> = OnceLock::new();
    let mask = *MASK.get_or_init(|| {
        std::env::var("LOSAT_PROBE_SHRINKMASK")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(u32::MAX)
    });
    if mask & bit != 0 {
        free(vector);
    } else {
        drop(vector);
    }
}

/// No NCBI counterpart: with LOSAT_TIMING, the `[TIMING] x_bn_shrink_free` line `main` prints at
/// exit (when a buffer was shrunk).
pub(crate) fn print_summary() {
    let shrunk = SHRUNK.load(Ordering::Relaxed);
    if crate::utils::xtiming::on() && shrunk > 0 {
        eprintln!("[TIMING] x_bn_shrink_free: buffers={shrunk}");
    }
}

/// A batch's BLASTNA query buffer, freed by `free` when the batch drops it.
pub(crate) struct Blastna(pub(crate) Vec<u8>);

impl std::ops::Deref for Blastna {
    type Target = Vec<u8>;

    fn deref(&self) -> &Vec<u8> {
        &self.0
    }
}

impl Drop for Blastna {
    fn drop(&mut self) {
        free_as(1, std::mem::take(&mut self.0));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_shrunk_vector_keeps_one_element_of_room() {
        // What `free` does before the drop, on vectors of every size class (heap and mapped).
        for &length in &[0usize, 1, 2, 100, 40_000, 1 << 20, 5 << 20] {
            let mut vector = vec![7u32; length];
            vector.clear();
            vector.shrink_to(1);
            assert!(vector.is_empty() && vector.capacity() <= length.max(1));
            free(vec![3u8; length]);
        }
    }
}
