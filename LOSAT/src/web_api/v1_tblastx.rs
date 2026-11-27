//! ABI v1's TBLASTX search (plan TD-1): the checks of the request stay frozen, NCBI's reader
//! reads the inputs (`v1_input`), and the records enter the engine's `run_local`. Moved from
//! `algorithm/tblastx/blast_engine/run_impl.rs` in session SFc, S10.

use anyhow::Result;

use super::v1_input::{read_queries, read_subjects};
use crate::algorithm::tblastx::blast_engine::run_local;
use crate::algorithm::tblastx::TblastxArgs;
use crate::api::local_blast::{OutputSink, ReportOutputs};

pub(super) fn run_web_pair(
    args: TblastxArgs,
    query_fasta: &str,
    subject_fasta: &str,
) -> Result<Vec<u8>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
    let subjects = read_subjects(subject_fasta.as_bytes(), "TBLASTX", false)?;
    let queries = read_queries(query_fasta.as_bytes(), "TBLASTX", false)?;
    // NCBI reference (598d8ae6): c++/include/algo/blast/api/local_blast.hpp:76-78
    // Web ABI v1 runs the same local search as the CLI and keeps the report in memory.
    let mut output = Vec::new();
    let outfmt = crate::utils::os_bytes::os_str_bytes(&args.outfmt).into_owned();
    let mut stderr = std::io::stderr();
    run_local(
        args,
        &queries,
        &crate::api::QueryEnd::Input,
        &subjects,
        &mut ReportOutputs::single(&outfmt, OutputSink::Writer(&mut output), &mut stderr),
    )?;
    Ok(output)
}
