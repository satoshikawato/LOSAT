//! `losat blastdbcmd -metadata`: the metadata of a database as NCBI's JSON
//! (`CSeqDB::GetDBMetaData` written by `CObjectOStreamJson`), the document makeblastdb writes
//! as `<db>.pjs` / `<db>.njs` (`writedb::metadata`).

use std::io::Write;

use crate::seqdb::date::iso8601_date_time;
use crate::seqdb::{SearchPath, SeqDb, SeqType};
use crate::writedb::metadata::{db_files_metadata, DbMetadata};

use super::args::Args;
use super::{database_error, refusal, Diag};

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1206-1222
// ```c++
//         else if (args["metadata"]) {
//         	x_InitBlastDB();
//             string output_prefix = args["metadata_output_prefix"]
//                     ? args["metadata_output_prefix"].AsString()
//                     : kEmptyStr;
//             if (!output_prefix.empty() && (output_prefix.back() != CFile::GetPathSeparator()))
//                 output_prefix += CFile::GetPathSeparator();
//         	CRef<CBlast_db_metadata> m = m_BlastDb->GetDBMetaData(output_prefix);
//             unique_ptr<CObjectOStreamJson> json_out(new CObjectOStreamJson(out, eNoOwnership));
//             json_out->SetDefaultStringEncoding(eEncoding_Ascii);
//             json_out->PreserveKeyNames();
//             CConstObjectInfo obj_info(m, m->GetTypeInfo());
//             json_out->WriteObject(obj_info);
//             json_out->Flush();
//             out.flush();
//             out << NcbiEndl;
//         }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1705-1747
// ```c++
// CRef<CBlast_db_metadata> CSeqDB::GetDBMetaData(string user_path)
// {
// 	CRef<CBlast_db_metadata> m (new CBlast_db_metadata());
// 	int num_seqs = 0;
// 	Uint8 total_length = 0;
//
// 	GetTotals(CSeqDB::eFilteredAll, &num_seqs, &total_length, true);
// 	vector<string> dblist;
//     NStr::Split(GetDBNameList(), " ", dblist, NStr::fSplit_Tokenize);
//     NON_CONST_ITERATE(vector<string>, itr, dblist) {
//         size_t off = (*itr).find_last_of(CFile::GetPathSeparator());
//         if (off != string::npos ) {
//         	(*itr).erase(0, off+1);
//         }
//     }
//
// 	string dbnames = NStr::Join(dblist, " ");
// 	m->SetDbname(dbnames);
//
//     m->SetDbtype(GetSequenceType() == CSeqDB::eProtein ? "Protein" : "Nucleotide" );
// 	m->SetDb_version(GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5?5:4);
// 	m->SetDescription(GetTitle());
// 	m->SetNumber_of_letters(total_length);
// 	m->SetNumber_of_sequences(num_seqs);
//
// 	CTimeFormat timeFmt = CTimeFormat::GetPredefined(CTimeFormat::eISO8601_DateTimeSec);
// 	CTime date(GetDate(), kBlastDbDateFormat);
// 	m->SetLast_updated(date.AsString(timeFmt));
//
// 	Int8 disk_bytes(0), cached_bytes(0);
// 	x_GetDBFilesMetaData(disk_bytes, cached_bytes, m->SetFiles(), user_path);
// 	m->SetBytes_total(disk_bytes);
// 	m->SetBytes_to_cache(cached_bytes);
//
// 	m->SetNumber_of_volumes(m_Impl->GetNumOfVols());
//
// 	if (GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5) {
// 		set<TTaxId> tax_ids;
// 		GetDBTaxIds(tax_ids);
// 		if((tax_ids.size() > 1) || ((tax_ids.size() == 1) && (0 != *tax_ids.begin()))){
// 			m->SetNumber_of_taxids(static_cast<int>(tax_ids.size()));
// 		}
// 	}
// ```
// The name is the `-db` value LOSAT opened as one database; a name with a double quote is
// refused (NCBI's name list would keep the quote in `dbname`). A date `CTime` does not read
// ends NCBI with an exception that names its source files (exit status 255, probes J189,
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
    out: &mut dyn Write,
) -> i32 {
    if name.contains('"') {
        return diag.refuse(&refusal(&format!(
            "-metadata of the database name {name} with a double quote"
        )));
    }
    let mut prefix = args
        .value("metadata_output_prefix")
        .unwrap_or("")
        .to_string();
    if !prefix.is_empty() && !prefix.ends_with(std::path::MAIN_SEPARATOR) {
        prefix.push(std::path::MAIN_SEPARATOR);
    }
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
        number_of_letters: db.total_length(),
        number_of_sequences: i64::from(db.num_seqs()),
        last_updated,
        number_of_volumes: db.volumes().len() as i32,
        number_of_taxids,
        bytes_total,
        bytes_to_cache,
        files,
    };
    // NCBI's stream ignores a failed write (no exception is enabled on it), so does this.
    let _ = out.write_all(&metadata.to_json());
    let _ = out.flush();
    0
}
