//! `losat blastdbcmd -metadata`: the metadata of a database as NCBI's JSON
//! (`CSeqDB::GetDBMetaData` written by `CObjectOStreamJson`), the document makeblastdb writes
//! as `<db>.pjs` / `<db>.njs` (`writedb::metadata`).

use crate::cli::NcbiStdio;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::OutputFile;
use crate::seqdb::date::iso8601_date_time;
use crate::seqdb::{SearchPath, SeqDb, SeqType};
use crate::writedb::cout::cout;
use crate::writedb::metadata::{db_files_metadata, metadata_user_path, write_metadata, DbMetadata};
use crate::writedb::ofstream::OfStream;

use super::{database_error, refusal, Diag};

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1206-1222
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1705-1747
// `GetDBNameList` is the `-db` value as given (`m_DBNames`, seqdbimpl.cpp:55): the names of a
// list go by spaces alone, and a double quote stays in the piece it is in (probes N121-N175 of
// line BDC). A date `CTime` does not read ends NCBI with an exception that names its source files (exit status 255, probes J189,
// J190 of line BDC), which LOSAT refuses. The JSON ends with the line feed of `NcbiEndl`
// (`writedb::metadata::write_metadata`, probes T13, T14, S39, J180-J199).
/// `-metadata` of the database `name` opened as `db`, its files named after
/// `-metadata_output_prefix`, to `out` (`-out`'s stream, nothing written to it yet); the exit
/// status.
pub(super) fn metadata(
    db: &SeqDb,
    name: &str,
    args: &Args,
    path: &SearchPath,
    diag: &mut Diag,
    out: OutputFile<'_>,
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
    // NCBI's serializer (`CObjectOStreamJson`) checks its writes and reports a failure where it
    // finds it: through its own buffer of 4096 bytes, on the `-out` file's `ofstream` or on
    // `cout` (`writedb::metadata::write_metadata`, `writedb::ofstream`, `writedb::cout`).
    // blastdbcmd 2.17.0 writing to a full device or under a file size limit gives these bytes
    // and exit statuses (line RG KC-11b-1). A closed pipe ends the program silently, where
    // SIGPIPE ends NCBI's (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10).
    let mut posts = Vec::new();
    let escaped = match out {
        OutputFile::File(file) => {
            let mut stream = OfStream::of_file(file);
            write_metadata(&mut stream, &metadata, super::TOOL.program, &mut posts)
        }
        OutputFile::Stdout(stdout) => {
            let mut stream = cout(stdout, NcbiStdio::cout());
            write_metadata(&mut stream, &metadata, super::TOOL.program, &mut posts)
        }
    };
    diag.write_text(&posts);
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
    match escaped {
        Some(message) => {
            diag.post(&format!("BLAST failed to write output: {message}"));
            6
        }
        None => 0,
    }
}
