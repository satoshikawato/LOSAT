//! EXPERIMENT (LOSAT_X_TBXRPTBATCH): the final order of TBLASTX's reports (`report_order` in
//! blast_engine/run_impl.rs) made for each query batch when the batch ends, instead of for every
//! batch at once after the last one.
//!
//! The reports are written once, after the last query batch, by the thread of the batch loop
//! (`BatchReports::write`). Before this switch that write began with the final order of every
//! batch's HSPs (`report::final_hit_order`, or `culled_hit_order` with `-culling_limit`): two
//! stable sorts per hit list, of HSPs of about 190 bytes each. Nothing else ran meanwhile: on the
//! genome 3 x 3 input at 16 threads (ws node 1, 2026-10-28, board 72960ee18) the report stage took
//! 0.187 s of 4.73 s, two thirds of it in that order.
//!
//! With the switch on, the batch loop puts each batch's HSPs in report order when it takes the
//! batch's result (`run_batch_loop`), and the write takes the list as it is. With the batch pipe
//! (TBXPIPE) the next batch is searched meanwhile on the other threads, so only the last batch's
//! order is left after the search. A `-db` search reads the titles of the batch's reported
//! sequences in report order (`database::read_titles`); it ordered a copy of the batch's HSPs for
//! that and now reads the batch's own ordered list.
//!
//! Why the reports are the same: the order of a query depends on its own HSPs only (its hit
//! lists, `-max_target_seqs`, `-max_hsps`, `-qcov_hsp_perc` and the culling of its subjects; see
//! `report_order`), and the batches hold disjoint, increasing ranges of queries. So the order of a
//! batch's queries made from that batch's HSPs is their order in the whole run, and the batches'
//! ordered lists one after the other, with the batch's query numbers shifted by its first query,
//! are the order of all of them. The `-db` path already relied on this for the titles.
//!
//! Off (`LOSAT_X_TBXRPTBATCH=0`): every batch's HSPs are kept as the search returns them and put
//! in order by the write, as before. Memory: the same list, held ordered instead of unordered;
//! the `-db` path no longer holds a copy of a batch's HSPs while it reads the titles.

use std::sync::OnceLock;

// No NCBI counterpart: chooses when the final order is made; it does not change any value NCBI
// computes. Read once per process.
/// LOSAT_X_TBXRPTBATCH (a strict-set switch: on unless `LOSAT_X_TBXRPTBATCH` is `0`;
/// `utils/xswitch.rs`).
pub(crate) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXRPTBATCH"))
}
