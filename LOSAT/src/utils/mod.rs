// No NCBI counterpart: the searches of a byte string that the ports share (cited in the module).
pub(crate) mod byte_search;
// No NCBI counterpart: decimal digits appended without `core::fmt` (cited in the module).
pub(crate) mod decimal;
pub mod dust;
pub mod entry_names;
// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp (`CFormatGuess`, cited in the module).
pub(crate) mod format_guess;
pub mod genetic_code;
pub mod glibc_fma_libm;
// NCBI reference (598d8ae6): c++/include/util/range_coll.hpp:331-334 (and the other calls of
// `lower_bound` named in the module): the position, asked in the order observed of
// `std::lower_bound`.
pub(crate) mod lower_bound;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:252-262 (`CRpsMmappedFile`,
// cited in the module).
pub(crate) mod mapped_file;
pub mod matrix;
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99, 615-621 (a file name's bytes in
// NCBI's messages, cited in the module).
pub(crate) mod os_bytes;
// No NCBI counterpart: the C library's text of an operating system error, which NCBI's messages
// end with, as NCBI BLAST+ 2.17.0 wrote it on Linux (cited in the module).
pub(crate) mod os_error_text;
pub(crate) mod process_id;
// No NCBI counterpart: whether two names are one file where there are no inode numbers (WASI).
pub(crate) mod same_file;
pub mod seg;
mod seg_lnfact;
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:665-674 (`-randomize`, cited
// at `writedb::raw_source::RawSeqDbSource::new`): `std::shuffle` with a `std::mt19937`, in the
// order observed of NCBI's Linux builds.
pub(crate) mod std_shuffle;
// No NCBI counterpart: test plumbing; the tests' temporary directory and process id
// (`std::env::temp_dir()` and `std::process::id()` panic on WASI), and the name NCBI's database
// reader gives such a directory (`test_dir::ncbi_dir`, cited there).
#[cfg(test)]
pub(crate) mod test_dir;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:177-188
pub mod threading;
pub mod winmask;
pub(crate) mod x87;
pub(crate) mod x_frame;
pub mod x_logclone;
pub(crate) mod x_mapped;
pub(crate) mod x_order;
pub mod x_seg_trim_memo;
pub(crate) mod xahead;
// No NCBI counterpart: item sets run by their owner and the pool's idle threads, no thread waiting
// in the pool (scheduling only, cited in the module).
pub mod xdiag;
pub(crate) mod xdrop_simd;
pub(crate) mod xenv;
pub(crate) mod xitems;
pub mod xreach;
pub mod xstats;
pub(crate) mod xswitch;
pub mod xtiming;
