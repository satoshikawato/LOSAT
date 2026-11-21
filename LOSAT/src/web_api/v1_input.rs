//! ABI v1's FASTA inputs of BLASTN, BLASTP and TBLASTX: NCBI's reader (`fasta_reader`) reads them
//! as it reads the files of the CLI and the registered inputs of ABI v2, and the warnings of
//! each record go with the record (`run_local` writes them). Plan TD-1, revised on 2026-11-05
//! (DW-48): the arguments, the formats and the errors of the requests stay frozen; the inputs
//! that ABI v1 accepts and what it reports about them are those of the CLI.

use anyhow::Result;

use crate::blastinput::app::engine_error;
use crate::blastinput::fasta_reader::{read_all, FastaInputSource, FastaRecord, ReaderConfig};
use crate::blastinput::input_files::is_blank;

/// The records of an input, as the CLI ends with a reader's error (`into_app_error`).
fn read_records(bytes: &[u8], config: ReaderConfig) -> Result<Vec<FastaRecord>> {
    read_all(
        &mut FastaInputSource::from_bytes(bytes, config),
        &mut |_| Ok(()),
    )
    .map_err(|error| error.into_app_error())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2554-2556
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:230-246
/// The subjects of `program`, which NCBI reads before the queries. An input without a record
/// is the engine's error (`run_local`, `CObjMgr_QueryFactory`).
pub(super) fn read_subjects(
    bytes: &[u8],
    program: &'static str,
    protein: bool,
) -> Result<Vec<FastaRecord>> {
    read_records(bytes, ReaderConfig::subject(program, protein, true))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:211-216
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:145-151
// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
/// The queries of `program`. A query of white space only has no record: NCBI's `Query is Empty!`
/// (`run_local`). A query of blank and comment lines only is not empty, and its first batch
/// (`eEOF`, no query) fails, as in `register` of ABI v2.
pub(super) fn read_queries(
    bytes: &[u8],
    program: &'static str,
    protein: bool,
) -> Result<Vec<FastaRecord>> {
    let records = read_records(bytes, ReaderConfig::query(program, protein, true))?;
    if records.is_empty() && !is_blank(bytes) {
        return Err(engine_error("Empty CBlastQueryVector"));
    }
    Ok(records)
}
