//! The lookup of an id in the tables of a database for `losat blastdbcmd`: NCBI's
//! `CSeqDBImpl` calls `AccessionToOids`, `AccessionsToOids` and `GiToOidwFilterCheck`, over the
//! LMDB files of a version 5 database (`CSeqDBLMDBSet`) or the ISAM files of the volumes of a
//! version 4 database (`CSeqDBVol::AccessionToOids`). The search programs look up the Seq-align
//! id of a reported sequence as NCBI's formatters find it again (`SeqidLookup`, `SeqidToOid`).

use std::sync::OnceLock;

use super::file::DbFile;
use super::filter::{entries_accessions_to_oids, LmdbEntry};
use super::isam::{seqdb_simplify_seqid, SeqDbIdType, SimplifiedAccession, StringIsam};
use super::isam_num::NumericIsam;
use super::lmdb_read::Env;
use super::volume::with_extension;
use super::{SeqDb, SeqDbError, SeqType};
use crate::seqid::{label, ParsedSeqId};

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2724-2728
/// The exception of a version 4 lookup whose number does not fit 32 bits.
pub const ID_OVERFLOWS: &str = "ID overflows range of specified type.";

/// The lookup tables of a database. A database restricted by a list (`SeqDb::restrict`) finds
/// only the OIDs the list includes, as NCBI's `x_CheckOrFindOID`; `losat blastdbcmd` restricts
/// none.
pub struct IdLookup<'a> {
    db: &'a SeqDb,
    /// The LMDB files of a version 5 database; none for version 4.
    entries: Vec<LmdbEntry>,
    isams: Isams<'a>,
    /// The first read of a candidate's header that `x_CheckVersions` makes (`take_utf8`).
    utf8: super::Utf8Reads,
}

/// The ISAMs of each volume (GIs; for version 4 also the strings and the PIGs), each opened at
/// the first lookup that reaches the volume.
pub(super) struct VolumeIsams {
    gi: Vec<OnceLock<Option<NumericIsam>>>,
    strings: Vec<OnceLock<Option<StringIsam>>>,
    pig: Vec<OnceLock<Option<NumericIsam>>>,
    /// Whether a string ISAM keeps its index file in memory (`StringIsam::keep_index`), as NCBI
    /// maps it: the search programs' lookups of the reported sequences (`SeqDb::seqid_lookup`).
    index_in_memory: bool,
}

impl VolumeIsams {
    fn new(volumes: usize, index_in_memory: bool) -> Self {
        Self {
            gi: cells(volumes),
            strings: cells(volumes),
            pig: cells(volumes),
            index_in_memory,
        }
    }
}

impl std::fmt::Debug for VolumeIsams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VolumeIsams")
            .field("volumes", &self.gi.len())
            .finish()
    }
}

/// The ISAMs of an `IdLookup`: its own (`losat blastdbcmd`) or the database's, kept for the
/// search (`SeqDb::seqid_lookup`).
enum Isams<'a> {
    Own(VolumeIsams),
    Shared(&'a VolumeIsams),
}

impl std::ops::Deref for Isams<'_> {
    type Target = VolumeIsams;

    fn deref(&self) -> &VolumeIsams {
        match self {
            Isams::Own(isams) => isams,
            Isams::Shared(isams) => isams,
        }
    }
}

impl SeqDb {
    /// The lookup tables of this database (for version 5, its LMDB files are read here).
    pub fn id_lookup(&self) -> Result<IdLookup<'_>, SeqDbError> {
        let entries = if self.is_version5() {
            self.lmdb_entries()?
        } else {
            Vec::new()
        };
        Ok(IdLookup {
            db: self,
            entries,
            isams: Isams::Own(VolumeIsams::new(self.volumes.len(), false)),
            utf8: Default::default(),
        })
    }
}

fn cells<T>(count: usize) -> Vec<OnceLock<T>> {
    (0..count).map(|_| OnceLock::new()).collect()
}

/// The value of `cell`, opened by `open` when the cell is empty.
fn opened<T>(
    cell: &OnceLock<Option<T>>,
    open: impl FnOnce() -> Result<Option<T>, SeqDbError>,
) -> Result<Option<&T>, SeqDbError> {
    if cell.get().is_none() {
        let _ = cell.set(open()?);
    }
    Ok(cell.get().and_then(Option::as_ref))
}

impl IdLookup<'_> {
    /// The first of the lookups' reads of a header (`x_CheckVersions` of an accession whose
    /// version the lookup stripped) that took a UTF8String for a VisibleString, which the caller
    /// takes after its lookup: NCBI posted its warning there (`seqdb::Utf8Reads`).
    pub fn take_utf8(&self) -> Option<super::Utf8Note> {
        self.utf8.take()
    }

    /// `x_CheckOrFindOID` of `oid` together with the test that `oid` is `oid0` once NCBI's OID list
    /// is set up (the search's lookups of its reported Seq-ids, after the search): an OID the
    /// database includes.
    fn visible(&self, oid: u32) -> bool {
        self.db.includes_oid(oid)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1237-1246
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:275-286
    /// Every OID of the key `acc` in a version 5 database: the LMDB files in the database's
    /// order, inside a file the order of LMDB's duplicates (`memcmp` of the 4 little-endian
    /// bytes). NCBI sets its OID list up first, so its check keeps the included OIDs.
    pub fn accession_to_oids(&self, acc: &[u8]) -> Result<Vec<u32>, SeqDbError> {
        self.db.build_oid_list()?;
        let mut oids = Vec::new();
        for entry in &self.entries {
            for oid in entry.oids(acc)? {
                if let Some(oid) = entry.database_oid(oid) {
                    if self.db.check_oid(oid)? {
                        oids.push(oid);
                    }
                }
            }
        }
        Ok(oids)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1347-1361
    /// One OID for each of `accs` in a version 5 database (`None`: not found), the first
    /// duplicate of the key in the last LMDB file that has it, kept when NCBI's check keeps it
    /// (`SeqDb::check_oid`: NCBI does not set its OID list up first, so the check of the first
    /// OID found can be made before). On an error `oids` holds what NCBI's vector holds at its
    /// exception: the answers of the LMDB files before the one that failed, unfiltered.
    pub fn accessions_to_oids(
        &self,
        accs: &[Vec<u8>],
        oids: &mut Vec<Option<u32>>,
    ) -> Result<(), SeqDbError> {
        entries_accessions_to_oids(&self.entries, accs, oids)?;
        for slot in oids.iter_mut() {
            if let Some(oid) = *slot {
                if !self.db.check_oid(oid)? {
                    *slot = None;
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1247-1279
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2791-2802
    /// Every OID of the user's id `acc` in a version 4 database, `simplified` its
    /// `SeqDB_SimplifyAccession` (which every volume computes alike): the volumes in order, in a
    /// volume the order of its ISAM (the lines of the string ISAM sort the OIDs of one key as
    /// decimal texts, `10` before `9`), an OID found before left out. `None` is NCBI's exception
    /// [`ID_OVERFLOWS`], after the lookup in the first volume.
    ///
    /// `gnl|BL_ORD_ID|n` is the one sequence `n` of the whole database (line BDC DECISIONS 12):
    /// NCBI adds `n` to the start of every volume. NCBI sets its OID list up first, so its
    /// check keeps the included OIDs.
    pub fn accession_to_oids_v4(
        &self,
        acc: &[u8],
        simplified: &SimplifiedAccession,
    ) -> Result<Option<Vec<u32>>, SeqDbError> {
        self.db.build_oid_list()?;
        if simplified.id_type == SeqDbIdType::Oid {
            let oid = u32::try_from(simplified.num_id).ok();
            return Ok(Some(match oid {
                Some(oid) if self.db.check_oid(oid)? => vec![oid],
                _ => Vec::new(),
            }));
        }
        let mut oids = Vec::new();
        for at in 0..self.db.volumes.len() {
            let vol_oids = self.volume_string_to_oids(at, acc, simplified)?;
            let ident = simplified.num_id;
            if !(ident == -1 || ident >> 32 == 0) {
                return Ok(None);
            }
            let vol_start = i64::from(self.db.volume_starts[at]);
            for oid in vol_oids {
                let Ok(oid) = u32::try_from(i64::from(oid) + vol_start) else {
                    continue;
                };
                if !self.db.check_oid(oid)? {
                    continue;
                }
                if !oids.contains(&oid) {
                    oids.push(oid);
                }
            }
        }
        Ok(Some(oids))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2638-2722
    /// `CSeqDBVol::x_StringToOids` of the volume `at` without its final checks (the size of the
    /// number, then `x_CheckVersions`, which [`StringIsam::string_id_to_oids`] runs): the OIDs
    /// local to the volume. A volume with a trace-id ISAM (`.nti`, which makeblastdb never
    /// writes) is refused when a `ti` id reaches it.
    fn volume_string_to_oids(
        &self,
        at: usize,
        acc: &[u8],
        simplified: &SimplifiedAccession,
    ) -> Result<Vec<i32>, SeqDbError> {
        let volume = &self.db.volumes[at];
        if volume.num_oids() == 0 {
            return Ok(Vec::new());
        }
        let seqtype = self.db.seqtype;
        let protein = seqtype == SeqType::Protein;
        // NCBI reads the index of a volume whose files fail `x_InitSearch` all the same
        // (`x_FindIndexBounds` in the constructor) and crashes (probe I of line BDC: an empty,
        // short or garbage index, another version or type, a data file of another length).
        let string_isam = || {
            opened(&self.isams.strings[at], || {
                let name = volume.base().to_string_lossy();
                let Some(isam) = StringIsam::open_existing(&name, protein)? else {
                    return Ok(None);
                };
                if !isam.is_initialized() {
                    return Err(SeqDbError::NcbiCrashes {
                        path: volume.base().to_path_buf(),
                        message: "a string ISAM whose header does not fit its files".to_string(),
                    });
                }
                if self.isams.index_in_memory {
                    return Ok(Some(isam.keep_index()?));
                }
                Ok(Some(isam))
            })
        };
        let local = |oids: Vec<u32>| oids.into_iter().map(|oid| oid as i32).collect();
        let found = |oid: Option<i32>| oid.into_iter().collect();
        Ok(match simplified.id_type {
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2768-2770
            // `x_CheckVersions` reads the header of each candidate, with the warning of a
            // UTF8String that its read takes for a VisibleString (`take_utf8`).
            SeqDbIdType::String => match string_isam()? {
                Some(isam) => local(isam.string_id_to_oids(acc, simplified, &mut |oid| {
                    let read = super::SeqLocation {
                        volume: at,
                        vol_oid: oid,
                        id: self.db.volume_starts[at] + oid,
                    };
                    // A read that fails may have posted the warning before (the note's path is
                    // then that of the string before the failure, `utf8_string_path`).
                    let set = volume
                        .defline_set(oid)
                        .inspect_err(|_| self.utf8.note(0, read))?;
                    if set.utf8 {
                        self.utf8.note(0, read);
                    }
                    Ok(set
                        .deflines
                        .into_iter()
                        .flat_map(|defline| defline.ids)
                        .collect())
                })?),
                None => Vec::new(),
            },
            SeqDbIdType::Pig => {
                let isam = opened(&self.isams.pig[at], || {
                    NumericIsam::open_pig_volume(volume.base(), seqtype)
                })?;
                match isam {
                    Some(isam) => found(isam.id_to_oid(i64::from(simplified.num_id as i32))?),
                    None => Vec::new(),
                }
            }
            SeqDbIdType::Gi => match self.gi_isam(at)? {
                Some(isam) => found(isam.id_to_oid(simplified.num_id)?),
                None => Vec::new(),
            },
            SeqDbIdType::Ti => {
                // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:214-223
                // The trace-id ISAM exists where `IndexExists` says so, as the string and PIG
                // ISAMs (`isam::index_exists`: two regular files).
                let index = with_extension(volume.base(), seqtype, "ti");
                if super::isam::index_exists(&index, &with_extension(volume.base(), seqtype, "td"))
                {
                    return Err(SeqDbError::Unsupported {
                        path: index,
                        message: "an ISAM of trace ids".to_string(),
                    });
                }
                match string_isam()? {
                    Some(isam) => local(isam.string_to_oids(acc, true, false)?.oids),
                    None => Vec::new(),
                }
            }
            SeqDbIdType::Oid => vec![simplified.num_id as i32],
        })
    }

    /// The GI ISAM of the volume `at` (`x_OpenGiFile`).
    fn gi_isam(&self, at: usize) -> Result<Option<&NumericIsam>, SeqDbError> {
        let volume = &self.db.volumes[at];
        if volume.num_oids() == 0 {
            return Ok(None);
        }
        opened(&self.isams.gi[at], || {
            NumericIsam::open_volume(volume.base(), self.db.seqtype)
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1190-1205
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1111-1122
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2413-2422
    // A volume's ISAM sets `oid` to -1 when it does not hold the PIG (`x_IdentToOid`); a
    // volume without PIG files (`x_OpenPigFile` opens none for an empty volume) leaves it
    // as it was, which in blastdbcmd is uninitialised (`oid` is declared without a value). No
    // list filters the OID.
    /// The OID of the PIG `pig`: `Some(Some(oid))` from the first volume whose PIG file has
    /// it, `Some(None)` when a volume has PIG files and none holds it (NCBI's -1), `None` when
    /// no volume has PIG files (NCBI's OID is then never set).
    pub fn pig_to_oid(&self, pig: i32) -> Result<Option<Option<u32>>, SeqDbError> {
        let mut searched = false;
        for (at, volume) in self.db.volumes.iter().enumerate() {
            if volume.num_oids() == 0 {
                continue;
            }
            let isam = opened(&self.isams.pig[at], || {
                NumericIsam::open_pig_volume(volume.base(), self.db.seqtype)
            })?;
            let Some(isam) = isam else {
                continue;
            };
            searched = true;
            if let Some(oid) = isam.id_to_oid(i64::from(pig))? {
                let oid = i64::from(oid) + i64::from(self.db.volume_starts[at]);
                return u32::try_from(oid).map(|oid| Some(Some(oid))).map_err(|_| {
                    SeqDbError::Format {
                        path: volume.base().to_path_buf(),
                        message: format!("the PIG file gives PIG {pig} the OID {oid}"),
                    }
                });
            }
        }
        Ok(searched.then_some(None))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1190-1206
    // NCBI's callers in blastdbcmd (`x_GetOids`, `x_ProcessBatchEntry_NoDup`) read `oid`, not
    // the answer: where the last volume holds the GI and the check does not keep its OID, `oid`
    // is where `CheckOrFindOID` left it (the next OID of the list, the range's start, or the
    // OID itself past the range's end), which blastdbcmd writes (cases A183-A190 of
    // tests/blastdbcmd_alias_keys.rs). A volume without the GI leaves -1 (`x_IdentToOid`).
    /// `GiToOidwFilterCheck` with `gi` and `oid`: NCBI's `oid` when it is not negative. The OID of
    /// the GI in the first volume whose file has it and whose OID NCBI's check keeps
    /// (`SeqDb::check_or_find_oid`; NCBI does not set its OID list up first), else, when the last
    /// volume has the GI, the OID that the check moved to.
    pub fn gi_to_oid(&self, gi: i64) -> Result<Option<u32>, SeqDbError> {
        let mut left = None;
        for at in 0..self.db.volumes.len() {
            left = None;
            let Some(isam) = self.gi_isam(at)? else {
                continue;
            };
            let Some(oid) = isam.id_to_oid(gi)? else {
                continue;
            };
            let Ok(found) = u32::try_from(i64::from(oid) + i64::from(self.db.volume_starts[at]))
            else {
                continue;
            };
            let mut oid = found;
            if self.db.check_or_find_oid(&mut oid)? && oid == found {
                return Ok(Some(oid));
            }
            left = Some(oid);
        }
        Ok(left)
    }
}

/// The LMDB files of a version 5 database, open for every [`SeqidLookup`] of a search
/// (`SeqDb::seqid_lookup`).
pub(super) struct LookupFiles {
    entries: Vec<LmdbEntry>,
    files: Vec<DbFile>,
}

impl std::fmt::Debug for LookupFiles {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LookupFiles")
            .field("files", &self.files.len())
            .finish()
    }
}

/// NCBI's lookup of a reported sequence's Seq-align id (`seqid_to_oid`), with the files of the
/// database kept open for the search (`SeqDb::seqid_lookup`): the LMDB files of version 5, here
/// as environments that keep the pages they read, and the ISAMs of the volumes, whose string
/// index stays in memory.
pub struct SeqidLookup<'a> {
    /// The ISAMs of the volumes: the strings of version 4, the GIs and the PIGs.
    volumes: IdLookup<'a>,
    /// Each LMDB file of the database with its environment (version 5), in the database's order.
    envs: Vec<(&'a LmdbEntry, Env<'a, DbFile>)>,
    /// What the lookups keep of each LMDB file of `envs` (`LmdbKept`).
    kept: std::cell::RefCell<Vec<LmdbKept>>,
    /// Whether the text-id lookups pass over an LMDB file without the accession table
    /// ([`SeqidLookup::in_indexed_files`]) instead of failing as NCBI's.
    skip_unindexed: bool,
}

/// What a `SeqidLookup` keeps of one LMDB file between its lookups.
#[derive(Default)]
struct LmdbKept {
    /// The `acc2oid` database once opened (`LmdbEntry::acc2oid_in`).
    acc2oid: Option<super::lmdb_read::Db>,
    /// Whether a key starts with each prefix that `case_folded_oids` asked about (the reported
    /// sequences' names share their first letters).
    prefixes: rustc_hash::FxHashMap<Vec<u8>, bool>,
}

impl LmdbKept {
    /// The `acc2oid` database of `entry`'s file opened as `env`, opened at the first lookup (an
    /// error is not kept: each lookup fails as NCBI's).
    fn acc2oid(
        &mut self,
        entry: &LmdbEntry,
        env: &Env<'_, DbFile>,
    ) -> Result<super::lmdb_read::Db, SeqDbError> {
        if let Some(db) = self.acc2oid {
            return Ok(db);
        }
        let db = entry.acc2oid_in(env)?;
        self.acc2oid = Some(db);
        Ok(db)
    }
}

impl SeqDb {
    /// The lookup of the Seq-align ids of the reported sequences of a search. The database's
    /// LMDB files (version 5) and ISAMs are opened at the first lookup of the search and stay
    /// open; each `SeqidLookup` reads the LMDB files through environments of its own.
    pub fn seqid_lookup(&self) -> Result<SeqidLookup<'_>, SeqDbError> {
        let mut envs = Vec::new();
        if self.is_version5() {
            let files = match self.lookup_files.get() {
                Some(files) => files,
                None => {
                    let entries = self.lmdb_entries()?;
                    let files = entries
                        .iter()
                        .map(LmdbEntry::open_file)
                        .collect::<Result<Vec<_>, _>>()?;
                    let _ = self.lookup_files.set(LookupFiles { entries, files });
                    self.lookup_files.get().expect("set above")
                }
            };
            for (entry, file) in files.entries.iter().zip(&files.files) {
                envs.push((entry, entry.open_env(file)?));
            }
        }
        let isams = self
            .lookup_isams
            .get_or_init(|| VolumeIsams::new(self.volumes.len(), true));
        Ok(SeqidLookup {
            volumes: IdLookup {
                db: self,
                entries: Vec::new(),
                isams: Isams::Shared(isams),
                utf8: Default::default(),
            },
            kept: std::cell::RefCell::new(envs.iter().map(|_| LmdbKept::default()).collect()),
            envs,
            skip_unindexed: false,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:190-198
    // `x_OpenPigFile` (seqdbvol.cpp:136-145, cited at `NumericIsam::open_pig_volume`) likewise;
    // the GI ISAM's failure to open is caught (`NumericIsam::open_volume`). The constructor of an
    // ISAM whose files exist (regular files, `CFile::Exists`) throws `Error: Could not open input
    // file` where one does not open (seqdbisam.cpp:1102-1171, cited at `StringIsam::open`), at
    // the first lookup that reaches the volume, which then fails on a file with NCBI's exception
    // (the prefetch's stop, `common::database::batch_prefetch`).
    /// Whether a lookup of a reported sequence's Seq-align id (`seqid_lookup`) may fail to open an
    /// ISAM of a volume with NCBI's exception: the open that the lookup makes of the string ISAM
    /// (`.psi`/`.psd`, `.nsi`/`.nsd`) or the PIG ISAM (`.ppi`/`.ppd`) of a volume with sequences,
    /// `isam::open_if_exists` (no ISAM where a file is no regular file, as NCBI's
    /// `CSeqDBIsam::IndexExists`; review 4 of line DBID loadord, H1), fails. Each file is opened
    /// and closed, never read (asked once for the run, `TaxdbRun::prefetch_can_stop`; review 3 of
    /// line DBID loadord, M1).
    pub fn isams_may_fail_to_open(&self) -> bool {
        let protein = self.seqtype == SeqType::Protein;
        let fails = |(index, data): (std::path::PathBuf, std::path::PathBuf)| {
            super::isam::open_if_exists(&index, &data).is_err()
        };
        self.volumes.iter().any(|volume| {
            if volume.num_oids() == 0 {
                return false;
            }
            StringIsam::file_names(&volume.base().to_string_lossy(), protein).map_or(true, fails)
                || fails(NumericIsam::pig_file_names(volume.base(), self.seqtype))
        })
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2679-2696
/// `IsStringId`: whether a version 5 database looks `id` up in its LMDB files (a `general` id
/// other than `PIG` falls through to the default).
fn is_string_id(id: &ParsedSeqId) -> bool {
    match id {
        ParsedSeqId::Gi(_) => false,
        ParsedSeqId::General(tag) => tag.db != b"PIG",
        _ => true,
    }
}

impl SeqidLookup<'_> {
    /// [`IdLookup::take_utf8`] of the lookups of the Seq-align ids.
    pub fn take_utf8(&self) -> Option<super::Utf8Note> {
        self.volumes.take_utf8()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:903-916
    /// `CSeqDB::SeqidToOid`: the OID that NCBI's lookup of the Seq-id `seqid` gives first, the
    /// one whose Bioseq NCBI's BLAST database data loader serves for it
    /// (`CBlastDbDataLoader::x_GetOid`, bdbloader.cpp:393; the search's own database, so its
    /// list's OIDs only), and its exception: the one lookup of NCBI's prefetch and formatters
    /// (`common::database::HeaderReader::found_instead`, `stop_lookup_unsure`). Without `multi`
    /// the ISAMs of the volumes return at the first OID found, so a volume after the first one
    /// that holds the id is never opened (an ISAM there that does not open stops none of them;
    /// review 3 of line DBID loadord, H3), where `seqid_to_oids` with `multi`, LOSAT's list of
    /// every holder, fails on it. `None`: not found (`seqid_to_oids`).
    pub fn seqid_to_oid(&self, seqid: &ParsedSeqId) -> Result<Option<u32>, SeqDbError> {
        Ok(self.seqid_to_oids(seqid, false)?.first().copied())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:178-190
    // makeblastdb writes the key of a Seq-id other than a text, `pdb`, `pir` or `prf` id also in
    // upper case, and for a `pdb` id its molecule and the molecule in upper case
    // (writedb_lmdb.cpp:136-150), where NCBI's lookup reads the key as given: `lcl|abc` does not
    // find `lcl|ABC`, which NCBI's scope takes for the same id (review 2 of line DBID N2b).
    /// The OIDs of a version 5 database under the upper-case key of `seqid` (for a `pdb` id its
    /// molecule in upper case), the search's list's OIDs only: every sequence that holds an id
    /// that differs from `seqid` in case alone. Empty for version 4, whose string ISAM keys are in
    /// lower case, and for the other ids (a text, `pir` or `prf` id: `case_variant_oids`).
    pub fn case_folded_oids(&self, seqid: &ParsedSeqId) -> Result<Vec<u32>, SeqDbError> {
        if !self.keys_in_lmdb(seqid)
            || seqid.textseq().is_some()
            || matches!(seqid, ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_))
        {
            return Ok(Vec::new());
        }
        let key = match seqid {
            ParsedSeqId::Pdb(pdb) => pdb.mol.to_ascii_uppercase(),
            _ => {
                let mut key = Vec::new();
                label::content(seqid, true, &mut key);
                key.make_ascii_uppercase();
                key
            }
        };
        self.lmdb_oids(&key)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:128-131
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:160-176
    // A text Seq-id's key has no upper-case copy: an accession's is in upper case
    // (`GetSeqIdString`, `seqid::label::content`), so the lookup finds every holder of it, but a
    // name alone and a `pir` or `prf` id's FASTA text (also its accession) are as written, and
    // NCBI's handles take `pir||pirname` and `pir||PIRNAME`, `gb||locusy` and `gb||LOCUSY` for one
    // (`CSeq_id_Textseq_Tree::x_FindStrInfo` and `x_Equals`, seq_id_tree.cpp:1064-1112;
    // `seqid::mask_handle::mask_handle_key`): LOSAT looks up every other key that differs in the
    // case of those letters (`tax_files::acc_case_variant_oids`; probe of line DBID pircase,
    // `blastdb-ids-20261021/R/pircase`). The key of a name alone is the same in every alternative
    // (`gb||locusy`, `emb||LOCUSY`, `sp||LOCUSY`), which the handles keep apart: the caller keeps
    // the sequences that hold an id of the same handle.
    /// The OIDs of a version 5 database under the other keys that differ from the key of `seqid`
    /// only in the case of its letters, the search's list's OIDs only, for a text id of a name
    /// alone and a `pir` or `prf` id (not the `pir|` or `prf|` before its FASTA text): the
    /// sequences that can hold an id that differs from `seqid` in case alone. Empty for version 4,
    /// whose string ISAM keys are in lower case, for a text id with an accession, whose key is in
    /// upper case, and for the other ids (`case_folded_oids`).
    pub fn case_variant_oids(&self, seqid: &ParsedSeqId) -> Result<Vec<u32>, SeqDbError> {
        if !self.keys_in_lmdb(seqid) {
            return Ok(Vec::new());
        }
        let mut key = Vec::new();
        let from = if matches!(seqid, ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_)) {
            label::fasta(seqid, &mut key);
            b"pir|".len()
        } else if seqid.textseq().is_some_and(|text| text.accession.is_none()) {
            label::content(seqid, true, &mut key);
            0
        } else {
            return Ok(Vec::new());
        };
        let mut oids = Vec::new();
        let mut kept = self.kept.borrow_mut();
        for ((entry, env), kept) in self.envs.iter().zip(kept.iter_mut()) {
            let Some(db) = self.acc2oid_of(kept, entry, env)? else {
                continue;
            };
            let found = entry.case_variant_oids_in(env, &db, &key, from, &mut kept.prefixes)?;
            oids.extend(
                found
                    .into_iter()
                    .filter_map(|oid| entry.database_oid(oid))
                    .filter(|&oid| self.volumes.visible(oid)),
            );
        }
        Ok(oids)
    }

    /// The OIDs of the key `key` in every LMDB file (`LmdbEntry::oids_of`), in the files' order,
    /// the search's list's OIDs only.
    fn lmdb_oids(&self, key: &[u8]) -> Result<Vec<u32>, SeqDbError> {
        let mut oids = Vec::new();
        let mut kept = self.kept.borrow_mut();
        for ((entry, env), kept) in self.envs.iter().zip(kept.iter_mut()) {
            let Some(db) = self.acc2oid_of(kept, entry, env)? else {
                continue;
            };
            let found = entry.oids_of(env, &db, key)?;
            oids.extend(
                found
                    .into_iter()
                    .filter_map(|oid| entry.database_oid(oid))
                    .filter(|&oid| self.volumes.visible(oid)),
            );
        }
        Ok(oids)
    }

    /// The `acc2oid` database of one LMDB file (`LmdbKept::acc2oid`): `None` for a file without
    /// the accession table where the lookup passes over such files ([`Self::in_indexed_files`]).
    fn acc2oid_of(
        &self,
        kept: &mut LmdbKept,
        entry: &LmdbEntry,
        env: &Env<'_, DbFile>,
    ) -> Result<Option<super::lmdb_read::Db>, SeqDbError> {
        match kept.acc2oid(entry, env) {
            Ok(db) => Ok(Some(db)),
            Err(error)
                if self.skip_unindexed
                    && error.ncbi_message().as_deref()
                        == Some(super::tax_files::NO_ACCESSION_INFO.as_bytes()) =>
            {
                Ok(None)
            }
            Err(error) => Err(error),
        }
    }

    /// Whether the accession table of one of the database's LMDB files cannot be had, so that a
    /// lookup of a text id fails on a file (`DB contains no accession info.` of a list that puts
    /// version 5 databases with and without the table side by side, NCBI's exception,
    /// `LookupFailure::ncbi_message`): the databases where the prefetch of a query batch can stop
    /// (`common::database::batch_prefetch`).
    pub fn lacks_an_accession_table(&self) -> bool {
        let mut kept = self.kept.borrow_mut();
        self.envs
            .iter()
            .zip(kept.iter_mut())
            .any(|((entry, env), kept)| kept.acc2oid(entry, env).is_err())
    }

    /// This lookup in the LMDB files that have the accession table only: where a list puts
    /// version 5 databases with and without it side by side, NCBI's lookup of a text id fails
    /// (`DB contains no accession info.`), but the sequences of the other files that hold such an
    /// id are those whose loads register it in NCBI's BLAST database data loader
    /// (`common::database::stop_lookup_unsure`). A file without the table holds no text id.
    pub fn in_indexed_files(mut self) -> Self {
        self.skip_unindexed = true;
        self
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1381-1459
    /// Whether [`SeqidLookup::seqid_to_oids`] looks `seqid` up by its key in the LMDB files (a
    /// text id of a version 5 database but `gnl|BL_ORD_ID|n`), not by `SeqDB_SimplifyAccession`.
    fn keys_in_lmdb(&self, seqid: &ParsedSeqId) -> bool {
        let bl_ord_id = matches!(seqid, ParsedSeqId::General(tag) if tag.db == b"BL_ORD_ID");
        self.volumes.db.is_version5() && !bl_ord_id && is_string_id(seqid)
    }

    /// The message of the exception that NCBI's lookup of `seqid` throws on the id itself, before
    /// it reads a file: `SeqDB_SimplifyAccession`'s (a string tag of `BL_ORD_ID` or `PIG`, a `ti`
    /// that is not a number, a `CSeqIdException`), which is no `CSeqDBException` (NCBI's
    /// prefetch catches it with its warning, `common::database::batch_prefetch`, where it rethrows
    /// the exception of a file). `None` where the lookup reads the files.
    pub fn throws_on_id(&self, seqid: &ParsedSeqId) -> Option<String> {
        if self.keys_in_lmdb(seqid) {
            return None;
        }
        seqdb_simplify_seqid(seqid)
            .err()
            .map(|error| error.to_string())
    }

    /// `CSeqDBImpl::SeqidToOids`: the OIDs that NCBI's lookup of the Seq-id `seqid` gives, in its
    /// order, the search's list's OIDs only; every one with `multi`, else the first.
    /// Version 5 looks a text id up by its key in every LMDB file (the files in order, in a file
    /// LMDB's order of the duplicates; every OID, as NCBI's), version 4 and the GIs and PIGs in the
    /// ISAMs of the volumes in order, and `gnl|BL_ORD_ID|n` is `n` from the start of each volume
    /// in turn. Empty: not found, also where NCBI's lookup throws on a number that does not fit 32
    /// bits. An error of a file is NCBI's exception of the lookup.
    pub fn seqid_to_oids(&self, seqid: &ParsedSeqId, multi: bool) -> Result<Vec<u32>, SeqDbError> {
        let db = self.volumes.db;
        if self.keys_in_lmdb(seqid) {
            let mut key = Vec::new();
            if matches!(seqid, ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_)) {
                label::fasta(seqid, &mut key);
            } else {
                label::content(seqid, true, &mut key);
            }
            return self.lmdb_oids(&key);
        }
        let simplified = seqdb_simplify_seqid(seqid)
            .map_err(|error| SeqDbError::NotFound(error.to_string().into_bytes()))?;
        let mut fasta = Vec::new();
        label::fasta(seqid, &mut fasta);
        let mut oids = Vec::new();
        for at in 0..db.volumes.len() {
            let vol_oids = self
                .volumes
                .volume_string_to_oids(at, &fasta, &simplified)?;
            let ident = simplified.num_id;
            if !(ident == -1 || ident >> 32 == 0) {
                return Ok(Vec::new());
            }
            let vol_start = i64::from(db.volume_starts[at]);
            for oid in vol_oids {
                if let Some(oid) = u32::try_from(i64::from(oid) + vol_start)
                    .ok()
                    .filter(|&oid| self.volumes.visible(oid))
                {
                    oids.push(oid);
                    if !multi {
                        return Ok(oids);
                    }
                }
            }
        }
        Ok(oids)
    }
}
