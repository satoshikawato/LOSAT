//! `-oid_masks` of makeblastdb: the OID mask file of each volume (`CWriteDB_OidList`) and the
//! test of a model sequence; the mask types and the file's name belong to `seqdb::oid_mask`,
//! which the reader of the alias key `OID_MASK_TYPE` shares, and the file's layout to
//! `seqdb::oid_list`, its reader.

use super::files::WriteDbFile;
use super::WriteDbError;
pub(crate) use crate::seqdb::oid_mask::EXCLUDE_MODEL;
use crate::seqdb::oid_mask::{exclude_model_extension, is_predicted};
use crate::seqid::accguide::SeqIdWarning;
use crate::seqid::ParsedSeqId;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:293-308
/// Whether `CWriteDB_Volume::WriteSequence` adds the OID of a sequence of the ids `ids` to the
/// list of models: every id that is not a GI is a predicted accession
/// (`IdentifyAccession()` has `fAcc_predicted` set). A sequence with no id but GIs is one, as in
/// makeblastdb 2.17.0 (line RG item D3b probe `p_v4`). Every id that is not a GI is identified,
/// also after one that is no model, and `warn` gets the accession guide's warnings of each
/// (review probe `sc`).
pub(crate) fn is_model(ids: &[ParsedSeqId], warn: &mut dyn FnMut(SeqIdWarning)) -> bool {
    let (models, accessions) = ids
        .iter()
        .filter(|id| !matches!(id, ParsedSeqId::Gi(_)))
        .fold((0usize, 0usize), |(models, accessions), id| {
            (models + usize::from(is_predicted(id, warn)), accessions + 1)
        });
    models == accessions
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:572-645
/// `CWriteDB_OidList` of `fExcludeModel`: the OIDs of a volume to mask, written at its close.
#[derive(Debug)]
pub(crate) struct OidList {
    file: WriteDbFile,
    oids: Vec<u32>,
}

impl OidList {
    pub(crate) fn exclude_model(dbname: &str, protein: bool, index: u32, max_fsize: u64) -> Self {
        Self {
            file: WriteDbFile::lazy(dbname, exclude_model_extension(protein), index, max_fsize),
            oids: Vec::new(),
        }
    }

    pub(crate) fn add_oid(&mut self, oid: u32) {
        self.oids.push(oid);
    }

    /// `Close`, given `total_oids`: the file of a volume of `total_oids` OIDs, none for 0.
    pub(crate) fn close(&mut self, total_oids: u32) -> Result<(), WriteDbError> {
        if total_oids != 0 && !self.file.created() {
            let mut map = vec![0xFFu8; ((total_oids - 1) / 8 + 1) as usize];
            for &oid in &self.oids {
                map[(oid >> 3) as usize] &= !(1u8 << (7 - (oid & 7)));
            }
            self.file.create_now();
            self.file.write_int4((total_oids - 1) as i32)?;
            self.file.write(&map)?;
            self.oids = Vec::new();
        }
        self.file.close();
        Ok(())
    }

    pub(crate) fn file_mut(&mut self) -> &mut WriteDbFile {
        &mut self.file
    }

    pub(crate) fn filename(&self) -> &str {
        self.file.filename()
    }
}

#[cfg(test)]
#[path = "oid_list_tests.rs"]
mod tests;
