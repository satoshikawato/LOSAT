//! A call boundary that keeps a closure's code out of its caller (code placement only).

// No NCBI counterpart: a call boundary; it computes no value NCBI computes.
/// Calls `f` in a function of its own that is never inlined: each closure gives its own instance,
/// so that the function that calls it keeps its own code (its registers, stack slots and branch
/// layout) and gains only the call. The callers name what they keep apart: the per-subject kernel
/// of BLASTP (another matrix than BLOSUM62, item SC2b), the scan of a TBLASTX chunk (the one-hit
/// loop of `-window_size 0`) and the unit search of RPS-BLAST and RPS-TBLASTN (the rest of a
/// query's search, series PR item PS10e).
#[inline(never)]
pub(crate) fn out_of_line<R>(f: impl FnOnce() -> R) -> R {
    f()
}
