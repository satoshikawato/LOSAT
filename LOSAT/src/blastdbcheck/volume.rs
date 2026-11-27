//! blastdbcheck's tests of a volume (`CTestActionList` of `CMetaDataTest`, `CStrideTest`,
//! `CSampleTest` and `CEndsTest`) and their check of one OID (`CTestAction::TestOID`).

use std::collections::BTreeSet;
use std::sync::Arc;

use crate::algorithm::common::database::TaxdbRun;
use crate::ncbi_app::Diag;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::matches_mask;
use crate::seqdb::lookup::SeqidLookup;
use crate::seqdb::rps::{read_cdd_freqs_header, read_cdd_rps_header, CddFreqsKind, RpsHeaderError};
use crate::seqdb::taxdb::LazyTaxDb;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqLocation, SeqType};
use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::DeflineSet;
use crate::utils::os_bytes::{path_bytes, path_from_bytes};
use crate::utils::random_gen::NcbiRandom;

use super::log::{Log, Stop, BRIEF, DETAILS, MINUTIAE};
use super::path::{entries_base, FileParts};
use super::work::{caught, unknown_text, Caught, Seen};
use super::PROGRAM;

/// The options that change what the tests check (`e_IsamLookup`, `e_Legacy`, `e_TaxIDSet`).
pub(super) struct Flags {
    /// The lookup of each Seq-id (not `-no_isam`).
    pub(super) isam: bool,
    /// The check for legacy files (`-legacy`).
    pub(super) legacy: bool,
    /// The check of a taxid for each sequence (`-must_have_taxids`).
    pub(super) taxids: bool,
}

/// State of `CCddDeltaHeadersTest` carried across all volumes of the run (never cleared between
/// volumes, matching NCBI BLAST+ 2.17.0).
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1185-1345
#[derive(Default)]
pub(super) struct CddDelta {
    magic_number: i32,
    num_profiles: i32,
    offsets: Vec<i32>,
    num_failures: u32,
}

/// A test of each volume, in the order the options add them.
pub(super) enum Test {
    /// `CMetaDataTest`, always first.
    MetaData,
    /// `CStrideTest` of every `N`-th OID (`-full`: 1, `-stride N`).
    Stride(i32),
    /// `CSampleTest` of `n` OIDs drawn by `random`, which every volume's draws go on with.
    Sample { n: i32, random: NcbiRandom },
    /// `CEndsTest` of `N` OIDs at each end.
    Ends(i32),
    /// `CCddDeltaHeadersTest` (`-cdd_delta`), checking `.rps`/`.freq`/`.wcounts`/`.obsr` headers.
    CddDelta(CddDelta),
}

impl Test {
    /// The test's name in its lines (`m_TestName`).
    fn name(&self) -> &'static [u8] {
        match self {
            Test::MetaData => b"MetaData",
            Test::Stride(_) => b"Stride",
            Test::Sample { .. } => b"Sample",
            Test::Ends(_) => b"EndPoints",
            Test::CddDelta(_) => b"CDD_Headers",
        }
    }
}

/// The tests of the run's volumes (`CTestActionList`) and the run's taxonomy database (NCBI's
/// static `CTaxDBFileInfo`, opened at the first lookup of any volume).
pub(super) struct Tests {
    list: Vec<Test>,
    flags: Flags,
    taxdb: Arc<LazyTaxDb>,
    path: SearchPath,
}

impl Tests {
    pub(super) fn new(list: Vec<Test>, flags: Flags, path: SearchPath) -> Self {
        let taxdb = Arc::new(LazyTaxDb::new(path.clone()));
        Tests {
            list,
            flags,
            taxdb,
            path,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:365-376
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:977-986
    /// The tests of the volume `name`, opened as `db` (`DoTests`): its faults. An exception that
    /// leaves a test (the OID list of `CheckOrFindOID`) ends the volume's tests with the driver's
    /// lines and one fault, whatever the tests before found.
    pub(super) fn run(
        &mut self,
        log: &mut Log<'_>,
        db: &SeqDb,
        name: &[u8],
        seen: &mut Seen,
        diag: &mut Diag<'_>,
    ) -> Result<u32, Stop> {
        let mut volume = Volume {
            db,
            name,
            flags: &self.flags,
            path: &self.path,
            taxdb: &self.taxdb,
            run: TaxdbRun::sequential_with(Arc::clone(&self.taxdb)),
            lookup: None,
            kept_set: DeflineSet::default(),
            diag,
        };
        let mut faults = 0;
        for test in &mut self.list {
            match volume.test(log, test, seen) {
                Ok(test_faults) => faults += test_faults,
                Err(Escape::Stop(stop)) => return Err(stop),
                Err(Escape::Thrown(error)) => {
                    let exception = caught(&error)?;
                    log.line(BRIEF, &[b"  [ERROR] caught exception in ", name])?;
                    log.what(b"", exception.what.as_deref(), &|| unknown_text(&error))?;
                    return Ok(1);
                }
            }
        }
        Ok(faults)
    }
}

/// What ends a test before it returns: the run's end, or an exception that the test does not
/// catch.
enum Escape {
    Stop(Stop),
    Thrown(SeqDbError),
}

impl From<Stop> for Escape {
    fn from(stop: Stop) -> Self {
        Escape::Stop(stop)
    }
}

/// The failure of the check of one OID: the step (`where`) and the exception.
struct Failure {
    step: &'static str,
    caught: Caught,
    /// LOSAT's refusal where the exception's text is written and LOSAT does not give it.
    unknown: String,
}

impl Failure {
    /// The failure of `step` with blastdbcheck's own exception (`std::runtime_error`) `text`.
    fn own(step: &'static str, text: impl Into<Vec<u8>>) -> Self {
        Failure {
            step,
            caught: Caught::own(text),
            unknown: String::new(),
        }
    }

    /// The failure of `step` where LOSAT's reader gives `error`; Err: LOSAT's refusal.
    fn of(step: &'static str, error: &SeqDbError) -> Result<Self, Stop> {
        Ok(Failure {
            step,
            caught: caught(error)?,
            unknown: unknown_text(error),
        })
    }
}

/// The site of `CSeqDBImpl::GetSeqIDs`'s throw for an OID past the volumes.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:813
const GET_SEQ_IDS_SITE: ThrowSite = ThrowSite {
    file: "src/objtools/blast/seqdb_reader/seqdbimpl.cpp",
    line: 813,
    kind: "CSeqDBException::eInvalid",
    function: Some("BLASTDB::ncbi::CSeqDBImpl::GetSeqIDs"),
};

/// The site of `CSeqDBImpl::GetTaxInfo`'s throw for a taxid that the taxonomy database does not
/// name.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1614-1622
const GET_TAX_INFO_SITE: ThrowSite = ThrowSite {
    file: "src/objtools/blast/seqdb_reader/seqdbimpl.cpp",
    line: 1620,
    kind: "CSeqDBException::eInvalid",
    function: Some("BLASTDB::ncbi::CSeqDBImpl::GetTaxInfo"),
};

/// The tests of one volume, opened as a database of its own (`CSeqDB` of its path), with that
/// database's caches: the taxonomy lookups of its volume (`run`) and its Seq-id lookup.
struct Volume<'a, 'd> {
    db: &'a SeqDb,
    /// The name the database was opened with (`GetDBNameList`).
    name: &'a [u8],
    flags: &'a Flags,
    path: &'a SearchPath,
    taxdb: &'a Arc<LazyTaxDb>,
    run: TaxdbRun,
    lookup: Option<SeqidLookup<'a>>,
    /// Kept header-set buffers reused across OIDs (`SeqDb::header_set_in`).
    kept_set: DeflineSet,
    diag: &'a mut Diag<'d>,
}

impl<'a> Volume<'a, '_> {
    /// The start of a test's line (`CTestAction::Log`): `  <db> / <test>: `.
    fn prefix(&self, test: &Test) -> Vec<u8> {
        [b"  ", self.name, b" / ", test.name(), b": "].concat()
    }

    /// Posts the `ERR_POST` texts of NCBI's reader.
    fn post(&mut self, texts: Vec<String>) {
        for text in texts {
            self.diag.err_post(&text);
        }
    }

    fn test(&mut self, log: &mut Log<'_>, test: &mut Test, seen: &mut Seen) -> Result<u32, Escape> {
        let prefix = self.prefix(test);
        match test {
            Test::MetaData => Ok(self.metadata(log, &prefix)?),
            &mut Test::Stride(n) => self.stride(log, &prefix, n, seen),
            Test::Sample { n, random } => self.sample(log, &prefix, *n, random, seen),
            &mut Test::Ends(n) => self.ends(log, &prefix, n, seen),
            Test::CddDelta(state) => Ok(self.cdd_delta(log, &prefix, state)?),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1185-1345
    /// `CCddDeltaHeadersTest::DoTest`: checks the `.rps` header against `.freq`, `.wcounts`, and
    /// `.obsr`, preserving `state` (offsets and failure count) across volumes as NCBI does.
    fn cdd_delta(
        &self,
        log: &mut Log<'_>,
        prefix: &[u8],
        state: &mut CddDelta,
    ) -> Result<u32, Stop> {
        let log_hdr_err = |log: &mut Log<'_>,
                           state: &mut CddDelta,
                           err: RpsHeaderError,
                           parts: &[&[u8]]|
         -> Result<(), Stop> {
            match err {
                RpsHeaderError::Ncbi(exception) => {
                    log.line(BRIEF, parts)?;
                    log.what(prefix, Some(&exception.what()), &|| String::new())?;
                    state.num_failures += 1;
                    Ok(())
                }
                RpsHeaderError::Refused(reason) => Err(Stop::Refused(reason)),
            }
        };
        match read_cdd_rps_header(self.name, self.path, PROGRAM) {
            Ok(header) => {
                state.magic_number = header.magic;
                state.num_profiles = header.num_profiles;
                state.offsets.extend_from_slice(&header.offsets);
            }
            Err(err) => log_hdr_err(
                log,
                state,
                err,
                &[
                    prefix,
                    b"  [ERROR] caught exception while reading rps header ",
                ],
            )?,
        }
        if state.num_failures != 0 {
            return Ok(state.num_failures);
        }
        for (kind, check_offsets) in [
            (CddFreqsKind::Freq, true),
            (CddFreqsKind::Wcounts, true),
            (CddFreqsKind::Obsr, false),
        ] {
            let ext = kind.ext().as_bytes();
            match read_cdd_freqs_header(self.name, self.path, kind, state.num_profiles, PROGRAM) {
                Ok(header) => {
                    if header.magic != state.magic_number {
                        log.line(
                            BRIEF,
                            &[
                                prefix,
                                b"  [ERROR] Bad magic number of the  ",
                                ext,
                                b" file",
                            ],
                        )?;
                        state.num_failures += 1;
                    }
                    if header.num_profiles != state.num_profiles {
                        let count = header.num_profiles.to_string();
                        log.line(
                            BRIEF,
                            &[
                                prefix,
                                b"  [ERROR] Wrong number of profiles in the  ",
                                ext,
                                b" file: ",
                                count.as_bytes(),
                            ],
                        )?;
                        state.num_failures += 1;
                    }
                    if check_offsets {
                        let count = state.num_profiles.max(0) as usize;
                        for i in 0..count {
                            if !matches!(
                                (header.offsets.get(i), state.offsets.get(i)),
                                (Some(a), Some(b)) if a == b
                            ) {
                                log.line(
                                    BRIEF,
                                    &[
                                        prefix,
                                        b"  [ERROR] Header offsets differ between .rps and ",
                                        ext,
                                        b" files",
                                    ],
                                )?;
                                state.num_failures += 1;
                            }
                        }
                    }
                }
                Err(err) => log_hdr_err(
                    log,
                    state,
                    err,
                    &[
                        prefix,
                        b"  [ERROR] caught exception while reading ",
                        ext,
                        b" freqs header ",
                    ],
                )?,
            }
        }
        Ok(state.num_failures)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:272-297
    /// `CheckOrFindOID` of `oid` (widened to `i64` per `Q-D2b-OV=A`: a negative one is moved to
    /// the first OID; one above `u32::MAX` is past every volume).
    fn check_or_find(&self, oid: &mut i64) -> Result<bool, Escape> {
        if *oid > i64::from(u32::MAX) {
            return Ok(false);
        }
        let mut next = u32::try_from(*oid).unwrap_or(0);
        let found = self
            .db
            .check_or_find_oid(&mut next)
            .map_err(Escape::Thrown)?;
        *oid = i64::from(next);
        Ok(found)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:392-470
    /// `CMetaDataTest`: its faults. Only the lookup of taxid 9606 can throw, before any line.
    fn metadata(&mut self, log: &mut Log<'_>, prefix: &[u8]) -> Result<u32, Stop> {
        let mut texts = Vec::new();
        let lookup = self.taxdb.lookup(9606, &mut texts);
        self.post(texts);
        let names = match lookup.and_then(|lookup| lookup.into_names()) {
            Ok(Some(names)) => names,
            Ok(None) => {
                let exception = NcbiException::new(GET_TAX_INFO_SITE, "Taxid 9606 not found");
                return metadata_exception(log, prefix, Caught::ncbi(&exception, true), &|| {
                    String::new()
                });
            }
            Err(error) => {
                let exception = caught(&error)?;
                return metadata_exception(log, prefix, exception, &|| unknown_text(&error));
            }
        };
        let db = self.db;
        let mut faults = 0;
        let mut fault = |log: &mut Log<'_>, parts: &[&[u8]]| -> Result<(), Stop> {
            faults += 1;
            let mut line = vec![prefix];
            line.extend_from_slice(parts);
            log.line(BRIEF, &line)
        };
        if db.date().is_empty() {
            fault(log, &[b"[ERROR] db has empty date string"])?;
        }
        let mut volumes: Vec<&std::path::Path> = db.volumes().iter().map(|v| v.base()).collect();
        volumes.sort_by(|a, b| SeqDb::compare_volume_paths(a, b));
        volumes.dedup();
        if volumes.is_empty() {
            fault(log, &[b"[ERROR] db has no volumes"])?;
        }
        for volume in volumes {
            for entry in volume_entries(&path_bytes(volume)) {
                let file = FileParts::of(&entry);
                if self.flags.legacy && matches_mask(file.ext, b"*tm") {
                    fault(log, &[b"[ERROR] legacy file ", &entry, b" exists."])?;
                }
                if file.file_length() == Some(0) {
                    fault(log, &[b"[ERROR] file ", &entry, b" has zero length."])?;
                }
            }
        }
        let oids = db.num_oids();
        let total = db.total_length();
        if i64::from(db.ncbi_num_seqs()) > i64::from(oids)
            || (total == 0) != (oids == 0)
            || (db.volume_length() == 0 && total != 0)
        {
            fault(log, &[b"[ERROR] sequence count/length mismatch"])?;
        }
        if names.scientific_name != b"Homo sapiens" {
            fault(
                log,
                &[
                    b"[ERROR] tax info looks wrong (",
                    &names.scientific_name,
                    b")",
                ],
            )?;
        }
        Ok(faults)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:592-603
    // Widen `oid + n` to `i64` so `oid` does not wrap past `INT_MAX`.
    /// `CStrideTest` of every `n`-th OID: 1 at the first OID that fails. With `n` of 0 or less
    /// NCBI's loop never ends after an OID that passes (`CheckOrFindOID` moves a negative OID to
    /// the first, which then passes as seen): LOSAT's refusal there.
    fn stride(
        &mut self,
        log: &mut Log<'_>,
        prefix: &[u8],
        n: i32,
        seen: &mut Seen,
    ) -> Result<u32, Escape> {
        let text = format!("<testing every {n}th OID>");
        log.line(MINUTIAE, &[prefix, text.as_bytes()])?;
        let mut oid = 0i64;
        while self.check_or_find(&mut oid)? {
            if !self.test_oid(log, prefix, seen, oid)? {
                return Ok(1);
            }
            if n <= 0 {
                return Err(Stop::Refused(format!(
                    "-stride {n} (NCBI's test of every {n}th OID never ends: it tests the first \
                     OID again and again), which is not supported by LOSAT's {PROGRAM}"
                ))
                .into());
            }
            oid += i64::from(n);
        }
        Ok(0)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:619-650
    /// `CSampleTest` of `n` OIDs drawn by `random`: 1 for a volume without OIDs, or at the first
    /// drawn OID that fails (in ascending order). NCBI draws until `n` OIDs pass
    /// `CheckOrFindOID`, which never ends where none does: LOSAT's refusal there.
    fn sample(
        &mut self,
        log: &mut Log<'_>,
        prefix: &[u8],
        n: i32,
        random: &mut NcbiRandom,
        seen: &mut Seen,
    ) -> Result<u32, Escape> {
        let max = i32::try_from(self.db.num_oids()).unwrap_or(i32::MAX);
        if max == 0 {
            log.line(BRIEF, &[prefix, b"  [ERROR] empty volume"])?;
            return Ok(1);
        }
        let mut first = 0i64;
        if n > 0 && !self.check_or_find(&mut first)? {
            return Err(Stop::Refused(format!(
                "a random sample of a database that includes none of its OIDs (NCBI draws OIDs \
                 without end), which is not supported by LOSAT's {PROGRAM}"
            ))
            .into());
        }
        let mut oids = BTreeSet::new();
        let mut drawn = 0;
        while drawn < n {
            let mut oid = i64::from(random.get_rand_range(0, (max - 1) as u32) as i32);
            if self.check_or_find(&mut oid)? {
                oids.insert(oid);
                drawn += 1;
            }
        }
        let text = format!(
            "<testing {n} randomly selected OIDs ({} unique)>",
            oids.len()
        );
        log.line(MINUTIAE, &[prefix, text.as_bytes()])?;
        for oid in oids {
            if !self.test_oid(log, prefix, seen, oid)? {
                return Ok(1);
            }
        }
        Ok(0)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:666-682
    // Approved exception (`PD-LOSAT-NCBI-DEFECTS`, Owner Q-D2b-OV=A): widen `count - n` and
    // `oid + 1` to `i64` so negative `-ends` near `INT_MIN` jump past the end of the volume as
    // nearby negative values (`-ends -1`, `-ends -2147483638`) do, rather than wrapping negative.
    /// `CEndsTest` of `n` OIDs at each end: 1 at the first OID that fails. From the `n`-th OID on
    /// the test jumps to the OID `n` before the end (the number of OIDs less `n`) where that is
    /// further, without `CheckOrFindOID`: with `n` of 0 or less past the last OID, which fails.
    fn ends(
        &mut self,
        log: &mut Log<'_>,
        prefix: &[u8],
        n: i32,
        seen: &mut Seen,
    ) -> Result<u32, Escape> {
        let text = format!("<testing {n} OIDs at each end>");
        log.line(MINUTIAE, &[prefix, text.as_bytes()])?;
        let count = i64::from(self.db.num_oids());
        let n = i64::from(n);
        let mut oid = 0i64;
        while self.check_or_find(&mut oid)? {
            if oid >= n {
                let end = count - n;
                if end > oid {
                    oid = end;
                }
            }
            if !self.test_oid(log, prefix, seen, oid)? {
                return Ok(1);
            }
            oid += 1;
        }
        Ok(0)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:474-582
    /// `TestOID`: whether `oid` passes (an OID that a test of the run checked before passes as
    /// seen), with its lines of Detailed and Minutiae.
    fn test_oid(
        &mut self,
        log: &mut Log<'_>,
        prefix: &[u8],
        seen: &mut Seen,
        oid: i64,
    ) -> Result<bool, Stop> {
        if !seen.insert(oid) {
            return Ok(true);
        }
        let failure = self.oid_checks(oid)?;
        if let Some(failure) = &failure {
            if log.shows(DETAILS) {
                let Some(what) = &failure.caught.what else {
                    return Err(Stop::Refused(failure.unknown.clone()));
                };
                let head = if failure.caught.seqdb {
                    format!("    Failed during {}, oid {oid}: ", failure.step)
                } else {
                    format!("    {} failed; oid={oid}: ", failure.step)
                };
                let text = [prefix, head.as_bytes(), what, b"\n"].concat();
                log.write(DETAILS, &text)?;
            }
        }
        let status = if failure.is_some() { "FAIL" } else { "PASS" };
        let text = format!("      Status for OID {oid}: {status}");
        log.line(MINUTIAE, &[prefix, text.as_bytes()])?;
        Ok(failure.is_none())
    }

    /// The steps of `TestOID` on `oid`: the failure of the first that fails.
    fn oid_checks(&mut self, oid: i64) -> Result<Option<Failure>, Stop> {
        let mut set = std::mem::take(&mut self.kept_set);
        let result = self.oid_checks_in(oid, &mut set);
        self.kept_set = set;
        result
    }

    /// [`Volume::oid_checks`] reusing `set` as its header buffer.
    fn oid_checks_in(&mut self, oid: i64, set: &mut DeflineSet) -> Result<Option<Failure>, Stop> {
        // Headers: `GetSeqIDs`, the Seq-ids of every definition line.
        let Some(at) = u32::try_from(oid)
            .ok()
            .and_then(|oid| self.db.locate_oid(oid).ok())
        else {
            let exception = NcbiException::new(GET_SEQ_IDS_SITE, "OID not found");
            return Ok(Some(Failure {
                step: "headers",
                caught: Caught::ncbi(&exception, true),
                unknown: String::new(),
            }));
        };
        let header = self.header_set(at, oid, set);
        // NCBI's `CObjectIStreamAsnBinary::FixVisibleChars` posts `Error: (803.7)` while
        // reading strings, both on valid headers and before a later BER error throws.
        let need_posts = match &header {
            Ok(true) => has_replaced_bytes(set),
            Ok(false) => false,
            Err(SeqDbError::Ncbi(_)) => true,
            Err(_) => false,
        };
        if need_posts {
            let posts = self.db.header_replaced_char_posts_at(at);
            self.diag.post_replaced(&posts);
        }
        let has_set = match header {
            Ok(has_set) => has_set,
            Err(error) => return Failure::of("headers", &error).map(Some),
        };
        // Sequence: `GetSequence`, which maps the sequence's bytes. For a protein with
        // non-ascending offsets, NCBI's `x_GetSequence` returns a negative length (which passes
        // the zero-length and null-pointer check here and then fails `GetBioseq` with no bioseq).
        let length = match self.sequence_len(at) {
            Ok(length) => length,
            Err(error) => return Failure::of("sequence", &error).map(Some),
        };
        if length == 0 {
            return Ok(Some(Failure::own("sequence", "sequence data is empty")));
        }
        // Bioseq: `GetBioseq`, none for no definition line or no letters (`length < 1`).
        if !has_set || set.deflines.is_empty() {
            return Ok(Some(Failure::own("bioseq", "no bioseq")));
        }
        if let Some(failure) = self.bioseq(at, oid as u32, set)? {
            return Ok(Some(failure));
        }
        if self.flags.taxids && !has_taxid(oid, set)? {
            return Ok(Some(Failure::own("taxid", "no taxid set")));
        }
        if self.flags.isam {
            return self.isam_lookups(oid, set);
        }
        Ok(None)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2086-2200,2250-2261
    /// `GetHdr` of the OID at `at` written over `set`: `false` for a stored header of no bytes, a
    /// set of no lines (`x_GetFilteredHeader`). Where the volume filters its sets NCBI's read of
    /// no set throws: LOSAT refuses it.
    fn header_set(
        &self,
        at: SeqLocation,
        oid: i64,
        set: &mut DeflineSet,
    ) -> Result<bool, SeqDbError> {
        if !self.db.header_set_in(at, set)? {
            if self.db.filters_header_sets(at) {
                return Err(SeqDbError::Unsupported {
                    path: path_from_bytes(self.name),
                    message: format!(
                        "OID {oid}: a stored header of no bytes in a volume whose sets NCBI filters"
                    ),
                });
            }
            return Ok(false);
        }
        Ok(true)
    }

    /// The length returned by `GetSequence` for the sequence at `at`: signed 32-bit for protein
    /// (`x_GetSequence`) and nucleotide (`x_GetSequence` / `GetSeqLengthExact`).
    fn sequence_len(&self, at: SeqLocation) -> Result<i64, SeqDbError> {
        match self.db.seqtype() {
            SeqType::Protein => self.db.protein_sequence_len_at(at).map(i64::from),
            SeqType::Nucleotide => self.db.nucleotide_length_i32_at(at).map(i64::from),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1306-1473
    /// The rest of `GetBioseq` of the OID at `at` whose definition lines are `set`: its length
    /// (none for no letters or negative signed 32-bit length), a nucleotide's letters with its
    /// ambiguities, and the taxonomy lookups of its sources, whose texts are posted.
    fn bioseq(
        &mut self,
        at: SeqLocation,
        oid: u32,
        set: &DeflineSet,
    ) -> Result<Option<Failure>, Stop> {
        let len = match self.db.seqtype() {
            SeqType::Protein => self.db.protein_length_i32_at(at).map(i64::from),
            SeqType::Nucleotide => self.db.nucleotide_length_i32_at(at).map(i64::from),
        };
        match len {
            Ok(n) if n < 1 => return Ok(Some(Failure::own("bioseq", "no bioseq"))),
            Ok(_) => {}
            Err(error) => return Failure::of("bioseq", &error).map(Some),
        }
        if self.db.seqtype() == SeqType::Nucleotide {
            let letters = self
                .db
                .read_nucleotide_search_sequence_at(at)
                .and_then(|block| block.nucleotide_seq_data(0).map(|_| ()));
            if let Err(error) = letters {
                return Failure::of("bioseq", &error).map(Some);
            }
        }
        match self.run.bioseq_texts(self.db, oid, set, None) {
            Ok(texts) => {
                self.post(texts);
                Ok(None)
            }
            Err(error) => Failure::of("bioseq", &error).map(Some),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:529-548
    /// The lookup of each Seq-id of `set` (`SeqidToOid`), which must give `oid`: the failure of the
    /// first that does not.
    fn isam_lookups(&mut self, oid: i64, set: &DeflineSet) -> Result<Option<Failure>, Stop> {
        const STEP: &str = "isam lookups";
        if self.lookup.is_none() {
            match self.db.seqid_lookup() {
                Ok(lookup) => self.lookup = Some(lookup),
                Err(error) => return Failure::of(STEP, &error).map(Some),
            }
        }
        let lookup = self.lookup.as_ref().expect("set above");
        for id in set.deflines.iter().flat_map(|defline| &defline.ids) {
            let found = match lookup.seqid_to_oid(id) {
                Ok(found) => found,
                Err(error) => return Failure::of(STEP, &error).map(Some),
            };
            let mut text = match found {
                None => b"seqid=".to_vec(),
                Some(found) if i64::from(found) != oid => {
                    format!("oid1={oid} oid2={found} seqid=").into_bytes()
                }
                Some(_) => continue,
            };
            crate::seqid::label::fasta(id, &mut text);
            return Ok(Some(Failure::own(STEP, text)));
        }
        Ok(None)
    }
}

/// The lines of an exception that `CMetaDataTest` caught: one fault.
fn metadata_exception(
    log: &mut Log<'_>,
    prefix: &[u8],
    exception: Caught,
    unknown: &dyn Fn() -> String,
) -> Result<u32, Stop> {
    log.line(BRIEF, &[prefix, b"  [ERROR] caught exception."])?;
    log.what(prefix, exception.what.as_deref(), unknown)?;
    Ok(1)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:425-453
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3958-4004
// `CDir::GetEntries` with the mask `<name>.???` (case sensitive) lists the directory in the order
// the system gives (`readdir`); a directory that cannot be read lists nothing.
/// The paths of the entries of the volume `volume`'s directory named `<volume's name>.???`.
fn volume_entries(volume: &[u8]) -> Vec<Vec<u8>> {
    let file = FileParts::of(volume);
    let name = &file.path[file.path.len() - file.base.len() - file.ext.len()..];
    let mask = [name, b".???"].concat();
    let base = entries_base(&file.dir);
    let Ok(entries) = std::fs::read_dir(path_from_bytes(base.as_slice())) else {
        return Vec::new();
    };
    entries
        .filter_map(Result::ok)
        .map(|entry| path_bytes(std::path::Path::new(&entry.file_name())).into_owned())
        .filter(|entry_name| matches_mask(entry_name, &mask))
        .map(|entry_name| [base.as_slice(), &entry_name].concat())
        .collect()
}

/// Whether any definition line in `set` holds bytes that NCBI's binary ASN.1 reader replaces with
/// `#` (and posts `Error: (803.7)` for).
fn has_replaced_bytes(set: &DeflineSet) -> bool {
    set.deflines
        .iter()
        .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:517-527
// NCBI reference (598d8ae6): c++/src/objects/seq/Bioseq.cpp:177-199
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1210-1231
// The Bioseq's taxid is that of its first source, and `GetBioseq` makes a source for each line of
// a positive taxid. Without one, `GetTaxidForSeqid` of the Bioseq's first Seq-id (the first line's
// first) takes the taxid of the first line whose first Seq-id of that type matches it
// (`CSeq_id::Match`, `compare`) and that has a taxid: the first line itself where it has one. A
// first line without Seq-ids makes NCBI read through a null reference.
/// Whether the sequence `oid` of the definition lines `set` has a taxid (`-must_have_taxids`).
fn has_taxid(oid: i64, set: &DeflineSet) -> Result<bool, Stop> {
    if set
        .deflines
        .iter()
        .any(|defline| defline.taxid.unwrap_or(0) > 0)
    {
        return Ok(true);
    }
    let Some(first_id) = set.deflines[0].ids.first() else {
        return Err(Stop::Refused(format!(
            "-must_have_taxids at OID {oid}, whose first definition line has no Seq-id (NCBI reads \
             a null reference there), which is not supported by LOSAT's {PROGRAM}"
        )));
    };
    for defline in &set.deflines {
        let same_type = defline
            .ids
            .iter()
            .find(|id| id.alternative() == first_id.alternative());
        if let (Some(id), Some(taxid)) = (same_type, defline.taxid) {
            if compare(id, first_id) == Sic::Yes {
                return Ok(taxid != 0);
            }
        }
    }
    Ok(false)
}

#[cfg(test)]
#[path = "volume_tests.rs"]
mod tests;
