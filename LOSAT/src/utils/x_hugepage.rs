//! EXPERIMENT (LOSAT_X_THP): ask the kernel for transparent huge pages on
//! large tables that are accessed at random (diagonal tables, lookup tables),
//! with `madvise(MADV_HUGEPAGE)`.  Placement only: no value changes.  Takes
//! effect when `/sys/kernel/mm/transparent_hugepage/enabled` is `madvise` or
//! `always` (Linux); elsewhere it is a no-op.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:143-150
//! ```c
//! /* Allocate the buffer to be used for diagonal array. */
//! diag_table->hit_level_array = (DiagStruct *)
//!     calloc(diag_table->diag_array_length, sizeof(DiagStruct));
//! if (word_params->options->window_size) {
//!     diag_table->hit_len_array = (Uint1 *)
//!          calloc(diag_table->diag_array_length, sizeof(Uint1));
//! }
//! ```
//! NCBI allocates the diagonal table with calloc, once per query length, and reads and writes it
//! at random. This module only advises the kernel how to back such a table in memory. It changes
//! no value and no access order, so every value NCBI computes stays the same (placement only).

use std::sync::OnceLock;

pub(crate) fn enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_THP").is_some())
}

/// Advises huge pages for the page-aligned interior of the buffer (at least
/// 4 MiB long); silently does nothing otherwise.
pub(crate) fn advise<T>(buffer: &mut [T]) {
    if !enabled() {
        return;
    }
    #[cfg(all(target_os = "linux", not(target_arch = "wasm32")))]
    {
        const PAGE: usize = 4096;
        const MADV_HUGEPAGE: i32 = 14;
        extern "C" {
            fn madvise(addr: *mut core::ffi::c_void, length: usize, advice: i32) -> i32;
        }
        let bytes = std::mem::size_of_val(buffer);
        if bytes < 4 << 20 {
            return;
        }
        let start = buffer.as_mut_ptr() as usize;
        let aligned_start = (start + PAGE - 1) & !(PAGE - 1);
        let aligned_end = (start + bytes) & !(PAGE - 1);
        if aligned_end > aligned_start {
            // SAFETY: the range lies inside the buffer and is page aligned; the
            // advice does not change the contents.
            unsafe {
                madvise(
                    aligned_start as *mut core::ffi::c_void,
                    aligned_end - aligned_start,
                    MADV_HUGEPAGE,
                );
            }
        }
    }
}
