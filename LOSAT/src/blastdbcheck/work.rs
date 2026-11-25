//! blastdbcheck's `Run` (`CBlastDbCheckApplication::Run`): the report's stream and header, the
//! tests that the arguments choose, and the walk of the volumes and alias files of the databases
//! that `-db` names (`CDbTest`) or that `-dir` holds (`CDirTest`), with its summary lines.

#[cfg(test)]
use std::cell::Cell;
use std::collections::{BTreeSet, HashSet};
use std::io::Write;

use crate::algorithm::common::database::losat_reader_error;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{Cin, OutputFile};
use crate::ncbi_app::Diag;
use crate::ncbi_exception::NcbiException;
use crate::seqdb::find::find_blast_dbs;
use crate::seqdb::taxdb::LazyTaxDb;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};
use crate::utils::random_gen::NcbiRandom;
use crate::writedb::time::LocalTime;

use super::alias::alias_file_test;
use super::log::{Log, Stop, BRIEF, DETAILS, SUMMARY};
use super::volume::{Flags, Test, Tests};
use super::PROGRAM;

/// The OIDs checked so far in the run (`TSeen`): an OID number met again, in any later volume, is
/// not checked again.
pub(super) type Seen = HashSet<i32>;

/// An exception that a `catch` of blastdbcheck takes: whether it is a `CSeqDBException`, and its
/// text (`what()`), `None` where LOSAT does not give NCBI's text of it.
pub(super) struct Caught {
    pub(super) seqdb: bool,
    pub(super) what: Option<Vec<u8>>,
}

impl Caught {
    /// An exception of the program's own (`std::runtime_error`) with the text `text`.
    pub(super) fn own(text: impl Into<Vec<u8>>) -> Self {
        Caught {
            seqdb: false,
            what: Some(text.into()),
        }
    }

    /// NCBI's exception `exception` (`CException`); `seqdb` whether it is a `CSeqDBException`.
    pub(super) fn ncbi(exception: &NcbiException, seqdb: bool) -> Self {
        Caught {
            seqdb,
            what: Some(exception.what()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:110-153
/// The exception that NCBI's reader throws where LOSAT's reader gives `error`, as a `catch` of
/// blastdbcheck takes it; `Err` where the error is LOSAT's (a database that LOSAT does not read,
/// a damaged file whose bytes NCBI reads past its end), whose outcome in NCBI's program LOSAT does
/// not know: LOSAT's refusal there. The text of a `CSeqDBException` is that of its place in NCBI's
/// reader (`SeqDbError::throw_site`), where LOSAT knows the place.
pub(super) fn caught(error: &SeqDbError) -> Result<Caught, Stop> {
    if let Some(exception) = error.ncbi_exception() {
        return Ok(Caught::ncbi(exception, false));
    }
    if !error.is_seqdb_exception() {
        return Err(Stop::Refused(losat_reader_error(error, PROGRAM)));
    }
    let what = match error {
        // NCBI's text of a failed integrity check names a line of its source; LOSAT's words do not.
        SeqDbError::Integrity { .. } => None,
        _ => error.ncbi_request_message().and_then(|message| {
            error
                .throw_site()
                .map(|site| NcbiException::new(site, message).what())
        }),
    };
    Ok(Caught { seqdb: true, what })
}

/// LOSAT's refusal where its report would write the text of `error`, which LOSAT does not give.
pub(super) fn unknown_text(error: &SeqDbError) -> String {
    format!(
        "the report of the error '{error}' at -verbosity 3 or more (LOSAT does not give NCBI's text \
         of it), which is not supported by LOSAT's {PROGRAM}"
    )
}

#[cfg(test)]
thread_local! {
    /// The clock's seconds that NCBI's sample test takes as its seed, as `CRandom::TValue`;
    /// the tests set it (`super::run_seeded`).
    pub(super) static SEED: Cell<Option<u32>> = const { Cell::new(None) };
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:612-617
/// The seed of the sample test: the current time in seconds, cut to 32 bits; `None` where the
/// platform has no clock.
fn time_seed() -> Option<u32> {
    #[cfg(test)]
    if let Some(seed) = SEED.with(Cell::get) {
        return Some(seed);
    }
    LocalTime::now().map(|time| time.unix_seconds as u32)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1376-1398
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-255
// The one exception that leaves the tests is the `CArgException` of `-logfile`'s
// `AsOutputFile`, which `CATCH_ALL` reports as `Command line argument error: <message>` with exit
// status 1 (probe m_bad_lf of row RG-D2b).
/// `CBlastDbCheckApplication::Run` on `args`; the exit status.
pub(super) fn run(
    args: &Args,
    settings: &ApplicationSettings,
    _stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    let mut file: OutputFile<'_>;
    let stream: &mut dyn Write = if args.given("logfile") {
        file = match args.output_file("logfile", out) {
            Ok(file) => file,
            Err(error) => {
                diag.post_bytes(
                    &[b"Command line argument error: ".as_slice(), error.message()].concat(),
                );
                return 1;
            }
        };
        if let Some(reason) = file.refusal() {
            return diag.refuse(&crate::seqmasks_io::run::unsupported(reason, PROGRAM));
        }
        &mut file
    } else {
        out
    };
    let level = args.integer("verbosity").unwrap_or(SUMMARY);
    let mut log = Log::new(stream, level);
    match check(args, settings, &mut log, diag) {
        Ok(status) => status,
        Err(Stop::Closed) => 0,
        Err(Stop::Refused(reason)) => diag.refuse(&reason),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:112-123
/// The name of a verbosity (`s_VerbosityString`).
fn verbosity_name(level: i32) -> &'static str {
    match level {
        0 => "Quiet",
        1 => "Brief",
        2 => "Summary",
        3 => "Detailed",
        4 => "Minutiae",
        _ => "",
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1386-1511
/// The work of `Run` after the report's stream is open: the header, the tests, the walk; the exit
/// status.
fn check(
    args: &Args,
    settings: &ApplicationSettings,
    log: &mut Log<'_>,
    diag: &mut Diag,
) -> Result<i32, Stop> {
    let db = args.bytes("db").unwrap_or_default();
    let dir = args.bytes("dir").unwrap_or_default();
    let target: Vec<u8> = match args.os_value("logfile") {
        Some(name) if args.given("logfile") => [
            b"file  (".as_slice(),
            name.as_encoded_bytes(),
            b")".as_slice(),
        ]
        .concat(),
        _ => b"<stdout>".to_vec(),
    };
    log.line(
        SUMMARY,
        &[
            b"Writing messages to ",
            &target,
            b" at verbosity (",
            verbosity_name(log.level()).as_bytes(),
            b")",
        ],
    )?;
    let dbtype = args.value("dbtype").unwrap_or("guess");
    if db.is_empty() == dir.is_empty() {
        log.line(BRIEF, &[b"error: Must specify exactly one of -dir or -db."])?;
        return Ok(1);
    }
    let flags = Flags {
        isam: !args.flag("no_isam"),
        legacy: args.flag("legacy"),
        taxids: args.flag("must_have_taxids"),
    };
    let enabled = |on: bool| if on { "EN" } else { "DIS" };
    for (what, on) in [
        ("ISAM", flags.isam),
        ("Legacy", flags.legacy),
        ("TaxID", flags.taxids),
    ] {
        log.line(
            SUMMARY,
            &[format!("{what} testing is {}ABLED.", enabled(on)).as_bytes()],
        )?;
    }
    let mut list = vec![Test::MetaData];
    if args.flag("full") {
        log.line(SUMMARY, &[b"Using `full' mode: every OID will be tested."])?;
        list.push(Test::Stride(1));
    }
    let sample = |n: i32| -> Result<Test, Stop> {
        let seed = time_seed().ok_or_else(|| {
            Stop::Refused(format!(
                "a random sample without a clock (its seed is the time), which is not supported \
                 by LOSAT's {PROGRAM}"
            ))
        })?;
        Ok(Test::Sample {
            n,
            random: NcbiRandom::new(seed),
        })
    };
    let mut chosen = args.flag("full");
    if let Some(stride) = args.integer("stride") {
        log.line(
            SUMMARY,
            &[format!("Testing every {stride}-th OID.").as_bytes()],
        )?;
        list.push(Test::Stride(stride));
        chosen = true;
    }
    if let Some(random) = args.integer("random") {
        log.line(
            SUMMARY,
            &[format!("Testing {random} randomly sampled OIDs.").as_bytes()],
        )?;
        list.push(sample(random)?);
        chosen = true;
    }
    if let Some(ends) = args.integer("ends") {
        log.line(
            SUMMARY,
            &[format!("Testing first {ends} and last {ends} OIDs.").as_bytes()],
        )?;
        list.push(Test::Ends(ends));
        chosen = true;
    }
    if !chosen {
        log.line(
            SUMMARY,
            &[b"By default, testing 200 randomly sampled OIDs."],
        )?;
        list.push(sample(200)?);
    }
    if args.flag("cdd_delta") {
        return Err(Stop::Refused(format!(
            "-cdd_delta (the check of a CDD database's headers for DELTA-BLAST), which is not \
             supported by LOSAT's {PROGRAM}"
        )));
    }
    log.line(SUMMARY, &[])?;
    let path = search_path(settings);
    let mut tests = Tests::new(list, flags, LazyTaxDb::new(path.clone()));
    let seqtype = match dbtype {
        "prot" => Some(SeqType::Protein),
        "nucl" => Some(SeqType::Nucleotide),
        _ => None,
    };
    let walk = Walk {
        path: &path,
        tests: &mut tests,
        diag,
    };
    let okay = if db.is_empty() {
        walk.directory(log, &dir, dbtype, args.flag("recursive"))?
    } else {
        walk.names(log, &db, seqtype)?
    };
    Ok(if okay { 0 } else { 1 })
}

/// The database search path of the run (`SeqDB_ResolveDbPath`'s: the current directory, `BLASTDB`
/// and the configuration's `[BLAST] BLASTDB`).
fn search_path(settings: &ApplicationSettings) -> SearchPath {
    SearchPath::from_environment(settings.blastdb.as_deref())
}

/// The walk of the run's volumes and alias files with its tests.
struct Walk<'r, 'd> {
    path: &'r SearchPath,
    tests: &'r mut Tests,
    diag: &'r mut Diag<'d>,
}

/// The volumes and alias files to test (`set<string>`: byte order, each once).
#[derive(Default)]
struct Found {
    volumes: BTreeSet<Vec<u8>>,
    aliases: BTreeSet<Vec<u8>>,
}

impl Walk<'_, '_> {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:927-1033
    // `NStr::Split` on single spaces (`fSplit_Tokenize`): runs of spaces part the names, quotes
    // are not read (probe r_space_db). With `guess`, NCBI's `FindVolumePaths` looks for each name
    // as a protein database and, after any exception, as a nucleotide one; each volume is opened
    // with NCBI's guess (`CSeqDB` of `eUnknown`).
    /// `CDbTest::Test` of the names of `db` of the type `seqtype` (`None`: `guess`): whether no
    /// test failed.
    fn names(
        mut self,
        log: &mut Log<'_>,
        db: &[u8],
        seqtype: Option<SeqType>,
    ) -> Result<bool, Stop> {
        let mut faults = 0u32;
        let mut found = Found::default();
        for name in db
            .split(|&byte| byte == b' ')
            .filter(|name| !name.is_empty())
        {
            faults += self.volume_list(log, name, seqtype, &mut found)?;
        }
        let total = found.volumes.len();
        log.line(SUMMARY, &[format!("Testing {total} volume(s).").as_bytes()])?;
        let mut seen = Seen::new();
        let mut passed = 0;
        for volume in &found.volumes {
            let volume_faults = self.volume(log, volume, seqtype, &mut seen)?;
            faults += volume_faults;
            passed += usize::from(volume_faults == 0);
        }
        Self::result(log, total, passed, "volume(s).")?;
        self.aliases_and_total(log, &found.aliases, faults)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1035-1178
    // The databases of the directory (`FindBlastDBs` with the aliases), each with its type: the
    // protein ones first, then the nucleotide ones, each volume opened as its type; `volumes.` in
    // the summary line where `CDbTest` says `volume(s).`.
    /// `CDirTest::Test` of the directory `dir`, `recursive`, of the databases of `dbtype`
    /// (`guess`: both types): whether no test failed.
    fn directory(
        mut self,
        log: &mut Log<'_>,
        dir: &[u8],
        dbtype: &str,
        recursive: bool,
    ) -> Result<bool, Stop> {
        log.line(SUMMARY, &[b"Finding database volumes."])?;
        let dbs = find_blast_dbs(dir, dbtype, recursive, false);
        let mut faults = 0u32;
        let mut protein = Found::default();
        let mut nucleotide = Found::default();
        for (seqtype, found) in [
            (SeqType::Protein, &mut protein),
            (SeqType::Nucleotide, &mut nucleotide),
        ] {
            for db in dbs.iter().filter(|db| db.seqtype == seqtype) {
                faults += self.volume_list(log, &db.quoted, Some(seqtype), found)?;
            }
        }
        let mut aliases = std::mem::take(&mut protein.aliases);
        aliases.append(&mut nucleotide.aliases);
        let total = protein.volumes.len() + nucleotide.volumes.len();
        log.line(SUMMARY, &[format!("Testing {total} volume(s).").as_bytes()])?;
        let mut seen = Seen::new();
        let mut passed = 0;
        for (seqtype, volumes) in [
            (SeqType::Protein, &protein.volumes),
            (SeqType::Nucleotide, &nucleotide.volumes),
        ] {
            for volume in volumes {
                let volume_faults = self.volume(log, volume, Some(seqtype), &mut seen)?;
                faults += volume_faults;
                passed += usize::from(volume_faults == 0);
            }
        }
        Self::result(log, total, passed, "volumes.")?;
        self.aliases_and_total(log, &aliases, faults)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:927-952
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1057-1078
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1038-1059
    // A catch of every exception: any exception of NCBI's skips the name; an error of LOSAT's is
    // its refusal.
    /// `x_GetVolumeList` of one name: its volumes and alias files into `found`; the number of
    /// faults (1 for a name whose files are not all found).
    fn volume_list(
        &mut self,
        log: &mut Log<'_>,
        name: &[u8],
        seqtype: Option<SeqType>,
        found: &mut Found,
    ) -> Result<u32, Stop> {
        let paths = |seqtype| SeqDb::volume_and_alias_paths(name, seqtype, self.path);
        let result = match seqtype {
            Some(seqtype) => paths(seqtype),
            None => match paths(SeqType::Protein) {
                Err(error) if ncbi_error(&error) => paths(SeqType::Nucleotide),
                result => result,
            },
        };
        match result {
            Ok((volumes, aliases)) => {
                found.volumes.extend(volumes);
                found.aliases.extend(aliases);
                Ok(0)
            }
            Err(error) if ncbi_error(&error) => {
                log.line(
                    BRIEF,
                    &[
                        b"  [ERROR] could not find all volume or alias files referenced in ",
                        name,
                        b", [skipped]",
                    ],
                )?;
                Ok(1)
            }
            Err(error) => Err(Stop::Refused(losat_reader_error(&error, PROGRAM))),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:968-986
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1104-1132
    /// The tests of the volume `volume` opened as `seqtype` (`None`: NCBI's guess); its faults.
    fn volume(
        &mut self,
        log: &mut Log<'_>,
        volume: &[u8],
        seqtype: Option<SeqType>,
        seen: &mut Seen,
    ) -> Result<u32, Stop> {
        log.line(DETAILS, &[b" ", volume])?;
        let types: &[SeqType] = match seqtype {
            Some(SeqType::Protein) => &[SeqType::Protein],
            Some(SeqType::Nucleotide) => &[SeqType::Nucleotide],
            None => &[SeqType::Protein, SeqType::Nucleotide],
        };
        let opened =
            SeqDb::open_guessing(types, |seqtype| SeqDb::open_in(volume, seqtype, self.path));
        match opened {
            Ok(db) => self.tests.run(log, &db, volume, seen, self.diag),
            Err(error) => {
                let exception = caught(&error)?;
                log.line(BRIEF, &[b"  [ERROR] caught exception in ", volume])?;
                log.what(b"", exception.what.as_deref(), &|| unknown_text(&error))?;
                Ok(1)
            }
        }
    }

    /// The summary line of `total` volumes or alias files of which `passed` passed; `unit` is
    /// what it counts.
    fn result(log: &mut Log<'_>, total: usize, passed: usize, unit: &str) -> Result<(), Stop> {
        let text = if total == passed {
            format!(" Result=SUCCESS. No errors reported for {total} {unit}")
        } else {
            format!(
                " Result=FAILURE. {} errors reported in {total} {unit}",
                total - passed
            )
        };
        log.line(SUMMARY, &[text.as_bytes()])
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1000-1033
    /// The test of each alias file of `aliases` (`CAliasTest`) and the run's last lines, with the
    /// faults `faults` of the names and volumes before: whether there was no fault.
    fn aliases_and_total(
        self,
        log: &mut Log<'_>,
        aliases: &BTreeSet<Vec<u8>>,
        mut faults: u32,
    ) -> Result<bool, Stop> {
        let total = aliases.len();
        log.line(SUMMARY, &[format!("Testing {total} alias(es).").as_bytes()])?;
        let mut passed = 0;
        for alias in aliases {
            log.line(DETAILS, &[b" ", alias])?;
            let alias_faults = alias_file_test(log, alias, self.path)?;
            faults += alias_faults;
            passed += usize::from(alias_faults == 0);
        }
        Self::result(log, total, passed, "alias(es).")?;
        if faults > 0 {
            log.line(BRIEF, &[])?;
            log.line(BRIEF, &[format!("Total errors: {faults}").as_bytes()])?;
        }
        Ok(faults == 0)
    }
}

/// Whether `error` is an exception of NCBI's reader (a `CSeqDBException` or another of NCBI's),
/// which a catch of every exception takes, rather than LOSAT's refusal.
fn ncbi_error(error: &SeqDbError) -> bool {
    error.is_seqdb_exception() || error.ncbi_exception().is_some()
}
