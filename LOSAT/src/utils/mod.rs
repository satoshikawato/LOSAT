pub mod dust;
pub mod genetic_code;
pub mod matrix;
pub mod seg;
mod seg_lnfact;
// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
// (*thread)->Run(); (*thread)->Join(&result);
pub mod threading;
