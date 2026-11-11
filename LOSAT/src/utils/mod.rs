pub mod dust;
pub mod entry_names;
pub mod genetic_code;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:252-262 (`CRpsMmappedFile`,
// quoted in the module).
pub(crate) mod mapped_file;
pub mod matrix;
pub mod seg;
mod seg_lnfact;
// libstdc++ (GCC 13): bits/stl_algo.h:1939-1951 (`std::sort`, quoted in the module): the order
// in which rpsbproc's sorts leave equal elements (RpsbProc/src/common/cdalignproc_base.cpp:200,
// ptrmap.hpp:141 of the 2025-05-13 package).
pub(crate) mod std_sort;
// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
// (*thread)->Run(); (*thread)->Join(&result);
pub mod threading;
pub mod winmask;
pub(crate) mod x87;
pub(crate) mod x_frame;
pub mod x_logclone;
pub(crate) mod x_mapped;
pub(crate) mod x_order;
pub mod x_seg_trim_memo;
pub(crate) mod xahead;
pub mod xdiag;
pub(crate) mod xdrop_simd;
pub(crate) mod xenv;
pub mod xreach;
pub mod xstats;
pub(crate) mod xswitch;
pub mod xtiming;
