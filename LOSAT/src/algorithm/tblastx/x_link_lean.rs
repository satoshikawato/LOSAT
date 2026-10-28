//! EXPERIMENT (LOSAT_X_TBXLINKLEAN): the working arrays of a TBLASTX linking call
//! (`link_group` in sum_stats_linking/linking_incr.rs) without the room they did not use.
//!
//! `link_group` holds about 200 bytes per HSP of its group while its rounds run, and several
//! groups link at once on the pool (on Sakai x MG1655 at 16 threads four groups of about 1.7 M
//! HSPs each). Two parts of that held more than the rounds read:
//!
//! 1. The maximum tree (`MaxTree`, 32 bytes per node) had a power of two of leaves, at least the
//!    group's HSPs, so a group of 2^k + 1 HSPs held twice the leaves it needed. With the switch it
//!    has one leaf per HSP after its inner nodes (`2 n` nodes, the leaf of HSP `i` at `n + i`). Each
//!    index in `2..2n` has its parent at `index / 2`, inside `1..n`, so node 1 is the root of every
//!    leaf for any `n`, and `set` and `root` read the same maxima: the keys are totally ordered
//!    (value, then list index), so the maximum does not depend on the tree's shape.
//! 2. The list of the index-1 HSPs (`elig1`) and the sorted distinct subject offsets (`so_vals`) are
//!    read only to build the sweep's orders and ranks before the rounds, and were held until the
//!    group's end. With the switch they are freed when the ranks are made.
//!
//! The rounds compute the same values in the same order: reports do not change.
//!
//! Off (`LOSAT_X_TBXLINKLEAN=0`): the tree has a power of two of leaves and the two lists live
//! until the group's end, as before.

use std::sync::OnceLock;

// No NCBI counterpart: the layout of LOSAT's own linking structures (NCBI's link_hsps.c keeps a
// linked list of HSPs and has no tree); it does not change any value NCBI computes. Read once per
// process.
/// LOSAT_X_TBXLINKLEAN (a strict-set switch: on unless `LOSAT_X_TBXLINKLEAN` is `0`;
/// `utils/xswitch.rs`).
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXLINKLEAN"))
}
