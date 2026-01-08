//! Filtering submodule for TBLASTX
//!
//! Contains HSP filtering functions.

mod purge_endpoints;

pub use purge_endpoints::purge_hsps_with_common_endpoints;
