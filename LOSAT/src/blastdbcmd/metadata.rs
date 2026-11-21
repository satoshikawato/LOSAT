//! `losat blastdbcmd -metadata`: the metadata of a database as NCBI's JSON
//! (`CSeqDB::GetDBMetaData` written by `CObjectOStreamJson`), the document makeblastdb writes
//! as `<db>.pjs` / `<db>.njs` (`writedb::metadata`).

use std::io::Write;

use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::Ostream;
use crate::seqdb::date::iso8601_date_time;
use crate::seqdb::{SearchPath, SeqDb, SeqType};
use crate::writedb::metadata::{db_files_metadata, metadata_user_path, DbMetadata};

use super::{database_error, refusal, Diag};

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1206-1222
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1705-1747
// `GetDBNameList` is the `-db` value as given (`m_DBNames`, seqdbimpl.cpp:55): the names of a
// list go by spaces alone, and a double quote stays in the piece it is in (probes N121-N175 of
// line BDC). A date `CTime` does not read ends NCBI with an exception that names its source files (exit status 255, probes J189,
// J190 of line BDC), which LOSAT refuses. The JSON ends with the line feed of `NcbiEndl`
// (`DbMetadata::to_json`, probes T13, T14, S39, J180-J199).
/// `-metadata` of the database `name` opened as `db`, its files named after
/// `-metadata_output_prefix`, to `out`; the exit status.
pub(super) fn metadata(
    db: &SeqDb,
    name: &str,
    args: &Args,
    path: &SearchPath,
    diag: &mut Diag,
    out: &mut Ostream<'_>,
) -> i32 {
    let prefix = metadata_user_path(args.value("metadata_output_prefix"));
    let dbname: Vec<&str> = name
        .split(' ')
        .filter(|piece| !piece.is_empty())
        .map(|piece| match piece.rfind(std::path::MAIN_SEPARATOR) {
            Some(off) => &piece[off + 1..],
            None => piece,
        })
        .collect();
    let last_updated = match iso8601_date_time(db.date()) {
        Ok(text) => text,
        Err(reason) => {
            return diag.refuse(&refusal(&format!(
                "{name}: NCBI does not read the date of the database ({reason}) and ends with \
                 an exception that names its source files (exit status 255)"
            )))
        }
    };
    let is_protein = db.seqtype() == SeqType::Protein;
    let (volumes, alias_files) = match SeqDb::volume_and_alias_paths(name, db.seqtype(), path) {
        Ok(paths) => paths,
        Err(error) => return database_error(diag, error),
    };
    let (bytes_total, bytes_to_cache, files) = db_files_metadata(
        &volumes,
        &alias_files,
        &db.lmdb_file_names(),
        is_protein,
        &prefix,
    );
    let number_of_taxids = if db.is_version5() {
        let tax_ids = match db.tax_id_lookup().and_then(|lookup| lookup.db_tax_ids()) {
            Ok(tax_ids) => tax_ids,
            Err(error) => return database_error(diag, error),
        };
        (tax_ids.len() > 1 || (tax_ids.len() == 1 && tax_ids.first() != Some(&0)))
            .then_some(tax_ids.len() as i32)
    } else {
        None
    };
    let metadata = DbMetadata {
        dbname: dbname.join(" ").into_bytes(),
        dbtype: if is_protein { "Protein" } else { "Nucleotide" },
        db_version: if db.is_version5() { 5 } else { 4 },
        description: db.title().to_vec(),
        number_of_letters: db.printed_total_length(),
        number_of_sequences: i64::from(db.printed_num_seqs()),
        last_updated,
        number_of_volumes: db.volumes().len() as i32,
        number_of_taxids,
        bytes_total,
        bytes_to_cache,
        files,
    };
    // NCBI's serial stream (`CObjectOStreamJson`) checks its writes: NCBI 2.17.0 with `-out
    // /dev/full` posts the failure of `json_out->Flush()` (`CObjectOStream: error at byte <the
    // JSON's length>`, three more lines of the stream's close, `BLAST failed to write output:
    // COStreamBuffer::Flush: failed`) and ends with exit status 6, or with other lines where a
    // JSON longer than the stream's buffer first fails (review 2 of line RG item K1, F2). LOSAT
    // refuses it, as dustmasker's serial writers do (RG DECISIONS K1-9). A closed pipe ends the
    // program silently, where SIGPIPE ends NCBI's (`Ostream::closed`).
    let _ = out.write_all(&metadata.to_json());
    let _ = out.flush();
    if out.bad() && !out.closed() {
        return diag.refuse(&format!(
            "a failed write of the -metadata output (NCBI's serial stream reports it) is not \
             supported by LOSAT's {}",
            super::PROGRAM
        ));
    }
    0
}
