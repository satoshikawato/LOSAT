//! `-oid_masks` of makeblastdb: the OID mask file of a volume (`CWriteDB_OidList`), `.pxm` or
//! `.nxm` for `EOidMaskType::fExcludeModel` (0x01), the one type makeblastdb 2.17.0 has.
//!
//! The file is the volume's last OID as a big-endian `Int4`, then one bit per OID, the first
//! OID of each byte in its most significant bit: 1 for an OID that stays, 0 for an OID whose
//! ids are all predicted RefSeq accessions (`XP_`, `XM_`, `XR_`, ...). A volume of no OID
//! writes no file. The reader of the alias key `OID_MASK_TYPE` reads these files.

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:147-150
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2707-2717

use super::files::WriteDbFile;
use super::WriteDbError;
use crate::seqid::accguide::{identify_seq_id, parse_flags, F_PREDICTED};
use crate::seqid::ParsedSeqId;

/// `EOidMaskType::fExcludeModel`.
pub(crate) const EXCLUDE_MODEL: u64 = 0x01;

/// `SeqDB_GetOidMaskFileExt` of `fExcludeModel`.
pub(crate) fn exclude_model_extension(protein: bool) -> &'static str {
    if protein {
        "pxm"
    } else {
        "nxm"
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:293-308
/// Whether `CWriteDB_Volume::WriteSequence` adds the OID of a sequence of the ids `ids` to the
/// list of models: every id that is not a GI is a predicted accession
/// (`IdentifyAccession() & fAcc_predicted`). A sequence with no id but GIs is one, as in
/// makeblastdb 2.17.0 (line RG item D3b probe `p_v4`).
pub(crate) fn is_model(ids: &[ParsedSeqId]) -> bool {
    ids.iter()
        .filter(|id| !matches!(id, ParsedSeqId::Gi(_)))
        .all(|id| {
            identify_seq_id(
                id,
                parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
                &mut |_| {},
            ) & F_PREDICTED
                != 0
        })
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

    /// `Close(total_oids)`: the file of a volume of `total_oids` OIDs, none for 0.
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
mod tests {
    use super::*;

    /// The bits of makeblastdb 2.17.0's `db.pxm` of 10 OIDs whose models are 0, 2, 4 and 9
    /// (probe `p_ps1`: `0000000957bf`), written to a file of its own.
    #[test]
    fn the_mask_file_is_makeblastdb_s() {
        let dir = std::env::temp_dir().join(format!("losat-oidlist-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let base = dir.join("db").to_string_lossy().into_owned();
        let mut list = OidList::exclude_model(&base, true, 0, 0);
        for oid in [0, 2, 4, 9] {
            list.add_oid(oid);
        }
        list.close(10).unwrap();
        let written = std::fs::read(list.filename()).unwrap();
        assert_eq!(written, [0, 0, 0, 9, 0x57, 0xBF]);
        let mut empty = OidList::exclude_model(&dir.join("e").to_string_lossy(), false, 0, 0);
        empty.close(0).unwrap();
        assert!(!std::path::Path::new(empty.filename()).exists());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
