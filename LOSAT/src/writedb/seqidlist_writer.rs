//! NCBI's writer of binary seqid lists (`WriteBlastSeqidlistFile`,
//! c++/src/objtools/blast/seqdb_writer/seqidlist_writer.cpp), which `blastdb_aliastool
//! -seqid_file_in` calls: the ids of a text list as a version 5 database looks them up
//! (`seqdb::seqidlist::process_for_v5`), kept to those of a database when one is named, in the
//! format that the readers of seqid lists read (`seqdb::seqidlist::write_binary`).

use super::time::LocalTime;
use crate::seqdb::alias::path_base_name;
use crate::seqdb::seqidlist::{self, SeqidListInfo};
use crate::seqdb::{SeqDb, SeqDbError};

/// A list that was not written.
#[derive(Debug)]
pub(crate) enum SeqidlistError {
    /// A `CSeqDBException` of the writer: its text, after `BLAST Database error:`.
    SeqDb(Vec<u8>),
    /// An id that does not parse: the message that NCBI posts (`ERR_POST` of the parser's
    /// exception), then the `CSeqDBException` text `Invalid seq id: <id>`.
    InvalidId { posted: Vec<u8>, message: Vec<u8> },
    /// No id is left after the gi ids are dropped: NCBI posts `Empty seqid list` and returns -1.
    Empty,
    /// The reader's error on the database's accessions.
    Reader(SeqDbError),
    /// LOSAT's error: the clock cannot be read for the list's date.
    Clock,
}

/// What the list says of the database `db` whose volume paths are `volumes` (`GetVolumeLength`,
/// the base names of `FindVolumePaths` with true joined by spaces, `GetDate`).
fn database_info(db: &SeqDb, volumes: &[String], info: &mut SeqidListInfo) {
    info.db_vol_length = db
        .volumes()
        .iter()
        .map(|volume| volume.index().volume_length)
        .sum();
    let names: Vec<&str> = volumes
        .iter()
        .map(|volume| path_base_name(volume))
        .collect();
    info.db_vol_names = names.join(" ").into_bytes();
    info.db_create_date = db.date().as_bytes().to_vec();
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/seqidlist_writer.cpp:67-187
// The database is checked first, then the ids are parsed (the warnings of a FASTA-style id of
// several parts are not shown at blastdb_aliastool's post level), and the empty list is posted
// before the database's accessions are looked up. The date is the clock's at the write
// (`kBlastDbDateFormat`).
/// The bytes of the binary list of the text lines `idlist` with `title` (`WriteBlastSeqidlistFile`
/// of the list, its title and the database `db`): the ids the database finds, with its length,
/// date and volumes, when `db` is given with its volume paths in NCBI's order
/// (`SeqDb::volume_and_alias_paths`).
pub(crate) fn seqidlist_bytes(
    idlist: &[Vec<u8>],
    title: &[u8],
    db: Option<(&SeqDb, &[String])>,
) -> Result<Vec<u8>, SeqidlistError> {
    let mut info = SeqidListInfo {
        title: title.to_vec(),
        ..SeqidListInfo::default()
    };
    if let Some((db, volumes)) = db {
        if !db.is_version5() {
            return Err(SeqidlistError::SeqDb(
                b"Seqidlist requires v5 database ".to_vec(),
            ));
        }
        database_info(db, volumes, &mut info);
    }
    let mut ids = seqidlist::process_for_v5(idlist, &mut |_| {}).map_err(|invalid| {
        SeqidlistError::InvalidId {
            posted: invalid.error.message.clone(),
            message: invalid.message(),
        }
    })?;
    if ids.is_empty() {
        return Err(SeqidlistError::Empty);
    }
    if let Some((db, _)) = db {
        let mut oids = Vec::new();
        db.id_lookup()
            .and_then(|lookup| lookup.accessions_to_oids(&ids, &mut oids))
            .map_err(SeqidlistError::Reader)?;
        ids = ids
            .into_iter()
            .zip(oids)
            .filter_map(|(id, oid)| oid.map(|_| id))
            .collect();
    }
    let now = LocalTime::now().ok_or(SeqidlistError::Clock)?;
    info.create_date = now.blast_db_date().into_bytes();
    Ok(seqidlist::write_binary(&info, &ids))
}
