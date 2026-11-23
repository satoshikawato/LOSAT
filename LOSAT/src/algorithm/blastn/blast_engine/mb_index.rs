//! `-use_index true` of a megablast search of a database: NCBI's loader of the MegaBLAST
//! database index (`CSetupFactory::InitializeMegablastDbIndex`, `DbIndexInit`), before the
//! formatter.
//!
//! NCBI looks for a new-style index first (the `-db` value: a super header `<volume>.shd` and
//! index volumes `<volume>.NN.idx` beside each database volume, `CIndexedDb_New`), then for an
//! old-style one (`-index_name`, or the `-db` value when the new style failed: index volumes
//! `<name>.NN.idx`, `CIndexedDb_Old`), and checks the word size against the first index volume
//! (`s_CheckMinWordSize`). Its error lines about the volumes it skips come first; when no index
//! loads, the first volume's header does not read, or the word size is below the index's, it
//! warns that the index will not be used and searches as without `-use_index` (it installs the
//! indexed search only while `-use_index` still holds), which LOSAT does with NCBI's bytes. When
//! an index loads, NCBI searches with it, as LOSAT does (`mb_index_seeds`, line RG items D5b-2
//! and D5b-3: also an index of only some database volumes, whose other volumes the search scans,
//! and the two-hit search; with `-db_hard_mask` the seeds placed into the hard ranges, item D5b-5,
//! the Owner's answer C to Q-HARDIDX), and LOSAT refuses the index files whose reading NCBI leaves
//! undefined or crashes on.

use std::path::Path;
use std::sync::Arc;

use anyhow::Result;

use super::database::BlastnDatabase;
use crate::algorithm::blastn::args::BlastnArgs;
use crate::dbindex::format::{volume_name, Header, SuperHeader, VolumeOids};
use crate::dbindex::search::{IndexedVolume, LoadedIndex};
use crate::dbindex::volume::IndexVolume;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::diag_log::{ERROR, WARNING};

/// What `-use_index true` gives the search (`megablast_index`).
#[derive(Debug)]
pub(super) enum MegablastIndex {
    /// NCBI searches without an index: it was not asked for, or none was used (after NCBI's
    /// warning).
    NotUsed,
    /// NCBI searches with an index that LOSAT does not search with, or reads one in a way that
    /// leaves its result undefined: LOSAT's refusal, the text after `-use_index true with `.
    Refused(String),
    /// NCBI searches with this index (`GetMBIndexLoaded`).
    Used(Arc<LoadedIndex>),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:2619-2632 (cited at
// `report::archive::options::algorithm_options`)
/// The program options of a search with an index: `CMbIndexArgs` set `ForceMbIndex` TRUE, which
/// stays when the index loads (the archive and the search strategy record FALSE for an index
/// that does not load, `algorithm_options`); set in the run's options handle, which the archive
/// and the search strategy share, before the batches.
pub(super) fn note_index_use(
    archive: Option<&crate::report::archive::ArchiveHandle>,
    strategy: Option<&crate::report::archive::strategy::StrategyExport>,
) -> Result<()> {
    let options = match (archive, strategy) {
        (Some(handle), _) => handle.options(),
        (None, Some(strategy)) => strategy.options(),
        (None, None) => return Ok(()),
    };
    options.set_search_sets(
        vec![crate::report::archive::options::Param {
            name: "ForceMbIndex",
            value: crate::report::archive::options::ParamValue::Boolean(true),
        }],
        Vec::new(),
        Vec::new(),
    )
}

/// LOSAT's refusal `-use_index true with {what}`.
pub(super) fn refusal(what: &str) -> anyhow::Error {
    anyhow::anyhow!("-use_index true with {what} is not supported by LOSAT's BLASTN")
}

impl MegablastIndex {
    /// The index that the search uses, or LOSAT's refusal of a search that NCBI runs with an
    /// index LOSAT does not search with (the loader's refusals). With `-db_hard_mask` the search
    /// places the index's seeds into the hard ranges (`mb_index_seeds::HardSeeds`; item D5b-5,
    /// the Owner's answer C to Q-HARDIDX).
    pub(super) fn check(&self) -> Result<Option<Arc<LoadedIndex>>> {
        match self {
            Self::NotUsed => Ok(None),
            Self::Refused(what) => Err(refusal(what)),
            Self::Used(index) => Ok(Some(Arc::clone(index))),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:216-221
// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:429-463
/// `-use_index true` of a megablast search of `-db` (`BlastnArgs::megablast_index_name`): NCBI's
/// loader, its error lines (`ERR_POST` at the error level) and its warning when it uses no index, each
/// one diagnostic (the archive saves each as one message, the warnings with code 1 and the error
/// lines with code 2). With a template (`-template_type` and `-template_length`, which require
/// each other) NCBI looks for no index.
pub(super) fn megablast_index(
    args: &BlastnArgs,
    database: &BlastnDatabase,
    diagnostics: &mut (dyn std::io::Write + Send),
) -> Result<MegablastIndex> {
    let Some((index_name, old_style)) = args.megablast_index_name() else {
        return Ok(MegablastIndex::NotUsed);
    };
    let (posts, result) = if args.template_length.is_some() {
        let text = "Database indexing is not available for discontiguous searches.";
        (Vec::new(), Init::NotUsed(text.as_bytes().to_vec()))
    } else {
        let word_size = crate::algorithm::blastn::coordination::determine_effective_word_size(args);
        let mut loader = Loader {
            database,
            posts: Vec::new(),
        };
        let result = loader.db_index_init(index_name, old_style, word_size);
        (loader.posts, result)
    };
    let index = match result {
        Init::NotUsed(text) => Err(text),
        Init::Used(loaded) => match loaded.into_index() {
            Ok(index) => Ok(index),
            Err(what) => return Ok(MegablastIndex::Refused(what)),
        },
        Init::Refused(what) => return Ok(MegablastIndex::Refused(what)),
    };
    for post in posts {
        let mut error = |text: &[u8]| {
            let line = [b"Error: [blastn] ".as_slice(), text, b"\n"].concat();
            write_post(args, diagnostics, ERROR, &line)
        };
        match post {
            Post::Line(text) => error(&text)?,
            Post::Gaps(from, to) => {
                for j in from..to {
                    error(format!("Index volume {j} not resolved.").as_bytes())?;
                }
            }
        }
    }
    let text = match index {
        Ok(index) => return Ok(MegablastIndex::Used(Arc::new(index))),
        Err(text) => text,
    };
    // The warning that `ERR_POST_EX` posts (code 1, subcode 1): `errstr`, then
    // " Database index will not be used.".
    let warning = [
        b"Warning: (1.1) [blastn] ".as_slice(),
        &text,
        b" Database index will not be used.\n",
    ]
    .concat();
    write_post(args, diagnostics, WARNING, &warning)?;
    Ok(MegablastIndex::NotUsed)
}

/// Writes NCBI's post `text` (of `severity`) in LOSAT's words, one message of the archive.
fn write_post(
    args: &BlastnArgs,
    diagnostics: &mut (dyn std::io::Write + Send),
    severity: i32,
    text: &[u8],
) -> Result<()> {
    let text = crate::wording::losat_bytes(text);
    match &args.archive {
        Some(handle) => handle.write_post(diagnostics, severity, &text)?,
        None => diagnostics.write_all(&text)?,
    }
    Ok(())
}

/// What `DbIndexInit` returns.
enum Init {
    /// NCBI's text of why no index is used (`errstr`, NCBI's words).
    NotUsed(Vec<u8>),
    /// NCBI searches with this index.
    Used(Loaded),
    /// NCBI's reading of the index is undefined, or it crashes: LOSAT's refusal.
    Refused(String),
}

/// NCBI's error lines (without `Error: [blastn] `).
enum Post {
    /// One line.
    Line(Vec<u8>),
    /// `Index volume j not resolved.` for each `j` of `from..to`, the old-style volumes missing
    /// between two found ones (as many as the range of `-index_name` has).
    Gaps(i64, i64),
}

/// The loader's state: the database and NCBI's error lines so far.
struct Loader<'a> {
    database: &'a BlastnDatabase,
    posts: Vec<Post>,
}

/// Why a loader did not load: NCBI's exception (`e.what()` goes into the result) or LOSAT's
/// refusal.
enum Failed {
    Ncbi(NcbiException),
    Refused(String),
}

/// What a loader found: the index volume whose header `s_CheckMinWordSize` reads, and the
/// volumes (`volumes_` of `CIndexedDb_New`: each index volume of an indexed database volume
/// with its number of sequences, or the path and the sequences of a database volume without an
/// index; the index volumes of `CIndexedDb_Old`, whose numbers of sequences their headers give).
struct Loaded {
    first: String,
    old_style: bool,
    volumes: Vec<Entry>,
}

/// One of `volumes_`: a path, whether it is an index volume, and its number of sequences (of
/// a new-style volume; an old-style one's is read from its header).
struct Entry {
    path: String,
    has_index: bool,
    n_oids: u32,
}

impl Loaded {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:812-857
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:1036-1066
    /// The index the search uses: its volumes with their first OIDs (`GetNextUnusedOID`; an
    /// old-style one's `seqmap_`, of the span from `StartSeq()` to `StopSeq()` of each header),
    /// each index volume read as the search reads it (`CDbIndex::Load`, which NCBI runs when the
    /// search reaches the volume), each database volume without an index kept for the scan (a
    /// mixed search). `Err`: LOSAT's refusal of a volume that does not read as an index.
    fn into_index(self) -> Result<LoadedIndex, String> {
        let mut volumes: Vec<IndexedVolume> = Vec::with_capacity(self.volumes.len());
        let mut next_oid: u32 = 0;
        for entry in self.volumes {
            let n_oids = if !entry.has_index {
                entry.n_oids
            } else {
                let volume = IndexVolume::open(Path::new(&entry.path)).map_err(|why| {
                    format!(
                        "a MegaBLAST index volume that LOSAT does not read ({}: {why})",
                        entry.path
                    )
                })?;
                if self.old_style {
                    volume.stop.wrapping_sub(volume.start)
                } else {
                    entry.n_oids
                }
            };
            volumes.push(IndexedVolume {
                path: entry.path,
                start_oid: next_oid,
                n_oids,
                has_index: entry.has_index,
            });
            next_oid = next_oid.wrapping_add(n_oids);
        }
        Ok(LoadedIndex::new(self.old_style, volumes))
    }
}

/// The place of NCBI's `CDbIndex_Exception` throws of `eBadOption` in `blast_dbindex.cpp`.
const fn dbindex_option(line: u32, function: &'static str) -> ThrowSite {
    ThrowSite {
        file: "src/algo/blast/api/blast_dbindex.cpp",
        line,
        kind: "CDbIndex_Exception::bad index creation option",
        function: Some(function),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:794-797
const NO_VOLUME_INDEX: ThrowSite =
    dbindex_option(795, "BLAST::ncbi::blast::CIndexedDb_New::CIndexedDb_New");
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:1023-1027
const NO_INDEX_FILE: ThrowSite =
    dbindex_option(1026, "BLAST::ncbi::blast::CIndexedDb_Old::CIndexedDb_Old");
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex.cpp:87-93
/// `CheckStream`'s throw for a stream that ended (`CIndexSuperHeader_Exception::eRead`).
const READ_FAILURE: ThrowSite = ThrowSite {
    file: "src/algo/blast/dbindex/dbindex.cpp",
    line: 92,
    kind: "CIndexSuperHeader_Exception::read failure",
    function: Some("ncbi::blastdbindex::{anonymous}::CheckStream"),
};

impl Loader<'_> {
    /// Keeps NCBI's error line `text`.
    fn error(&mut self, text: String) {
        self.posts.push(Post::Line(text.into_bytes()));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:1122-1171
    /// `DbIndexInit` of `indexname`, `old_style`, `partial` and `word_size`: the new style (unless
    /// `-index_name` asks for the old one), then the old style after a failure, whose
    /// exceptions' texts join with a line feed; a loaded index gets the word size check, whose
    /// text replaces any earlier one.
    fn db_index_init(&mut self, indexname: &str, old_style: bool, word_size: usize) -> Init {
        let mut result = Vec::new();
        if !old_style {
            match self.new_style(indexname) {
                Ok(loaded) => return self.check_min_word_size(loaded, word_size),
                Err(Failed::Ncbi(e)) => result = e.ncbi_what(),
                Err(Failed::Refused(what)) => return Init::Refused(what),
            }
        }
        match self.old_style(indexname) {
            Ok(loaded) => self.check_min_word_size(loaded, word_size),
            Err(Failed::Ncbi(e)) => {
                result.push(b'\n');
                result.extend_from_slice(&e.ncbi_what());
                // `what()` ends with a line feed, which NCBI's post shows before " Database
                // index will not be used.".
                Init::NotUsed(result)
            }
            Err(Failed::Refused(what)) => Init::Refused(what),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:753-803
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:631-651
    /// `CIndexedDb_New` of `indexname` and `partial`: the database of the `-db` value (its names split
    /// at single spaces and joined again by `SeqDB_CombineAndQuote`, which drops only leading
    /// spaces that `CSeqDB` ignores anyway), its volumes in `FindVolumePaths` order, each with
    /// `AddIndexInfo`; an exception when no volume has an index, or the database reader's when
    /// the value does not open as it is (white space other than spaces at its ends, which the
    /// search's own opening drops).
    fn new_style(&mut self, indexname: &str) -> Result<Loaded, Failed> {
        let volumes = self
            .database
            .index_volumes(indexname)
            .map_err(|error| match error {
                Some(exception) => Failed::Ncbi(exception),
                None => Failed::Refused(format!(
                    "a -db value that NCBI's index loader opens apart from the search ('{indexname}')"
                )),
            })?;
        // `volumes_`: each index volume of an indexed database volume, or the database
        // volume's own path for one without an index.
        let mut entries: Vec<Entry> = Vec::new();
        for (vol_name, dbnseq) in &volumes {
            self.add_index_info(vol_name, *dbnseq, &mut entries)?;
        }
        if !entries.iter().any(|entry| entry.has_index) {
            return Err(Failed::Ncbi(NcbiException::new(
                NO_VOLUME_INDEX,
                "no database volume has an index",
            )));
        }
        // `MinIndexWordSize` reads `volumes_[0].name`: an index volume, or the path of a
        // first database volume without an index, which does not read as an index header
        // ("Failed to read ...": no index is used, Q-FIRSTVOL).
        Ok(Loaded {
            first: entries[0].path.clone(),
            old_style: false,
            volumes: entries,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:664-750
    /// `AddIndexInfo` (with `partial`) of the database volume `vol_name` of `dbnseq`
    /// sequences: its index volumes into `entries` when its index is used, else the volume
    /// itself after NCBI's error lines (a missing or unreadable super header is posted at the
    /// level `Info`, which is not shown). A volume without sequences whose super header lists
    /// no index volume adds nothing.
    fn add_index_info(
        &mut self,
        vol_name: &str,
        dbnseq: u32,
        entries: &mut Vec<Entry>,
    ) -> Result<(), Failed> {
        let resolved = self.resolve_volume_index(vol_name, dbnseq)?;
        match resolved {
            Some(names) => entries.extend(names.into_iter().map(|(path, n_oids)| Entry {
                path,
                has_index: true,
                n_oids,
            })),
            None => entries.push(Entry {
                path: vol_name.to_string(),
                has_index: false,
                n_oids: dbnseq,
            }),
        }
        Ok(())
    }

    /// `AddIndexInfo`'s checks of the index of `vol_name`: its index volumes with their numbers
    /// of sequences, or `None` (not resolved) after the error line of the check that failed.
    fn resolve_volume_index(
        &mut self,
        vol_name: &str,
        dbnseq: u32,
    ) -> Result<Option<Vec<(String, u32)>>, Failed> {
        let Some(shdr) = SuperHeader::load(Path::new(&format!("{vol_name}.shd"))) else {
            return Ok(None);
        };
        if shdr.num_seq != dbnseq {
            self.error(format!(
                "numbers of OIDs reported by the database and by the index do not match. \
                     Index for volume {vol_name} will not be used"
            ));
            return Ok(None);
        }
        let mut names = Vec::new();
        let mut total_idxvol_oids: u64 = 0;
        for i in 0..shdr.num_vol {
            let Some(name) = self.database.resolve_db_path(&volume_name(vol_name, i)) else {
                self.error(format!(
                    "index volume  not resolved; index will not be used for {vol_name}"
                ));
                return Ok(None);
            };
            let idxvol_oids = match crate::dbindex::format::volume_oids(Path::new(&name)) {
                VolumeOids::Count(count) => count,
                VolumeOids::ReadFailure(place) => {
                    return Err(Failed::Ncbi(NcbiException::new(
                        READ_FAILURE,
                        format!("[{name}] {place}"),
                    )))
                }
                VolumeOids::Undefined => {
                    return Err(Failed::Refused(format!(
                        "a MegaBLAST index volume whose header NCBI does not read whole ({name})"
                    )))
                }
            };
            if idxvol_oids == 0 {
                self.error(format!(
                    "index volume {name} reports no sequences; index will not be used for \
                         {vol_name}"
                ));
                return Ok(None);
            }
            names.push((name, idxvol_oids));
            total_idxvol_oids += u64::from(idxvol_oids);
        }
        if u64::from(dbnseq) != total_idxvol_oids {
            self.error(format!(
                "total of oids reported by index volumes ({total_idxvol_oids}) does not match \
                     the number of oids reported by the superheader ({dbnseq}); index will not be \
                     used for {vol_name}"
            ));
            return Ok(None);
        }
        Ok(Some(names))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:944-1028
    /// `CIndexedDb_Old` of `indexnames`: the names (separated by runs of spaces) in
    /// `SeqDB_CompareVolume` order, each `base[,ignored[,first[,last]]]` (volumes 0 to 99 by
    /// default, `atoi` of the numbers), the volumes `<base>.NN.idx` found on the search path
    /// (`SeqDB_ResolveDbPath`), an error line for each volume missing between two found ones;
    /// an exception when none is found.
    fn old_style(&mut self, indexnames: &str) -> Result<Loaded, Failed> {
        let mut index_names: Vec<String> = Vec::new();
        if !indexnames.is_empty() {
            let mut dbnames = split_old_style_names(indexnames);
            dbnames.sort_by(|a, b| {
                crate::seqdb::SeqDb::compare_volume_paths(Path::new(a), Path::new(b))
            });
            for indexname in &dbnames {
                let (index_base, start_vol, stop_vol) = parse_old_style_name(indexname);
                if start_vol > stop_vol {
                    continue;
                }
                // NCBI's loop counts `i` (a `long`) up by one from `start_vol` while `i` cast to
                // `unsigned long` is at most `stop_vol`: a last volume of `ULONG_MAX`, or one
                // past `LONG_MAX` after a first one below it, never ends the loop (the `long`
                // overflows).
                if stop_vol == u64::MAX
                    || (stop_vol > i64::MAX as u64 && start_vol <= i64::MAX as u64)
                {
                    return Err(Failed::Refused(format!(
                        "an -index_name range of volumes that NCBI's loader does not finish ('{indexname}')"
                    )));
                }
                let mut last_i: i64 = -1;
                let mut i = start_vol as i64;
                while (i as u64) <= stop_vol {
                    let number = if i >= 0 {
                        format!("{i:02}")
                    } else {
                        i.to_string()
                    };
                    if let Some(name) = self
                        .database
                        .resolve_db_path(&format!("{index_base}.{number}.idx"))
                    {
                        if last_i + 1 < i {
                            self.posts.push(Post::Gaps(last_i + 1, i));
                        }
                        index_names.push(name);
                        last_i = i;
                    }
                    i += 1;
                }
            }
        }
        match index_names.first() {
            Some(first) => Ok(Loaded {
                first: first.clone(),
                old_style: true,
                volumes: index_names
                    .into_iter()
                    .map(|path| Entry {
                        path,
                        has_index: true,
                        n_oids: 0,
                    })
                    .collect(),
            }),
            None => Err(Failed::Ncbi(NcbiException::new(
                NO_INDEX_FILE,
                format!("no index file specified or index '{indexnames}*' not found."),
            ))),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:1098-1120
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:453-457
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:200-204
    /// `s_CheckMinWordSize` of `word_size` on a loaded index: the minimum of the index volume
    /// `loaded.first`. A word size below it resets the index, and a header that does not read
    /// gives NCBI's "Failed to read" text with the index still loaded; after either text NCBI
    /// clears `-use_index` and so never installs the indexed search. An empty volume crashes
    /// NCBI's reading.
    fn check_min_word_size(&self, loaded: Loaded, word_size: usize) -> Init {
        let min_ws = match min_index_word_size(Path::new(&loaded.first)) {
            MinWordSize::Read(min_ws) => min_ws,
            MinWordSize::Unreadable => {
                return Init::NotUsed(b"Failed to read index MegaBLAST db min word size.".to_vec())
            }
            MinWordSize::Empty => {
                return Init::Refused(format!(
                    "an empty MegaBLAST index volume ({})",
                    loaded.first
                ))
            }
        };
        if (word_size as i64) < i64::from(min_ws) {
            return Init::NotUsed(
                format!("MegaBLAST database index requires word size greater than {min_ws}.")
                    .into_bytes(),
            );
        }
        Init::Used(loaded)
    }
}

/// What `s_GetMinimumSupportedWordSizeByIndex` reads of a file.
enum MinWordSize {
    /// The minimum word size of the header.
    Read(i32),
    /// NCBI's mapping or header check throws.
    Unreadable,
    /// An empty regular file, which NCBI maps as no memory and reads through a null pointer.
    Empty,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:906-942
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5747-5771
/// `s_GetMinimumSupportedWordSizeByIndex` of `fname`: the header of the file `fname` as `CMemoryFile`
/// maps it. The mapping throws for a path that is not a regular file (after links) or does not
/// open, maps an empty file as nothing, and reads zeros past the end of a file shorter than the
/// header; the header check throws for a word size field over 15. The path is examined before it
/// is opened, as NCBI does (opening a FIFO would wait).
fn min_index_word_size(fname: &Path) -> MinWordSize {
    let Ok(metadata) = std::fs::metadata(fname) else {
        return MinWordSize::Unreadable;
    };
    if !metadata.is_file() {
        return MinWordSize::Unreadable;
    }
    if metadata.len() == 0 {
        return MinWordSize::Empty;
    }
    let Ok(file) = std::fs::File::open(fname) else {
        return MinWordSize::Unreadable;
    };
    let mut bytes = Vec::with_capacity(crate::dbindex::format::HEADER_SIZE);
    let read = std::io::Read::read_to_end(
        &mut std::io::Read::take(file, crate::dbindex::format::HEADER_SIZE as u64),
        &mut bytes,
    );
    match read.ok().and_then(|_| Header::read(&bytes)) {
        Some(header) => MinWordSize::Read(header.min_word_size()),
        None => MinWordSize::Unreadable,
    }
}

/// `CIndexedDb_Old`'s names: the text up to each space, then from the next character that is
/// not a space (`find_first_of(" ")`, `find_first_not_of(" ")`): a leading space gives an
/// empty first name, trailing spaces none.
fn split_old_style_names(indexnames: &str) -> Vec<String> {
    let bytes = indexnames.as_bytes();
    let mut names = Vec::new();
    let mut start = Some(0usize);
    while let Some(at) = start {
        let end = bytes[at..].iter().position(|&b| b == b' ').map(|o| at + o);
        names.push(indexnames[at..end.unwrap_or(bytes.len())].to_string());
        start = end.and_then(|end| {
            bytes[end..]
                .iter()
                .position(|&b| b != b' ')
                .map(|o| end + o)
        });
    }
    names
}

/// `CIndexedDb_Old`'s reading of one name `base[,ignored[,first[,last]]]`: the base and the
/// first and last volumes (0 and 99 unless given; `atoi` of the text, kept as `unsigned long`).
fn parse_old_style_name(indexname: &str) -> (String, u64, u64) {
    let bytes = indexname.as_bytes();
    let len = bytes.len();
    let find = |from: usize| -> Option<usize> {
        bytes
            .get(from..)
            .and_then(|rest| rest.iter().position(|&b| b == b',').map(|o| from + o))
    };
    // `substr(start, end - start)` with `end` npos: to the end.
    let substr = |from: usize, end: Option<usize>| -> &[u8] {
        &bytes[from.min(len)..end.unwrap_or(len).max(from.min(len))]
    };
    let atoi = |text: &[u8]| -> u64 { crate::ncbi_str::c_atoi(text) as i64 as u64 };
    let mut start_vol: u64 = 0;
    let mut stop_vol: u64 = 99;
    let end = find(0);
    let index_base = String::from_utf8_lossy(substr(0, end)).into_owned();
    // `start = end + 1`: npos + 1 is 0.
    let mut start = end.map_or(0, |end| end + 1);
    if start < len && end.is_some() {
        let end = find(start);
        start = end.map_or(0, |end| end + 1);
        if start < len && end.is_some() {
            let end = find(start);
            let start_vol_str = substr(start, end);
            if !start_vol_str.is_empty() {
                start_vol = atoi(start_vol_str);
            }
            start = end.map_or(0, |end| end + 1);
            if start < len && end.is_some() {
                let end = find(start);
                let stop_vol_str = substr(start, end);
                if !stop_vol_str.is_empty() {
                    stop_vol = atoi(stop_vol_str);
                }
            }
        }
    }
    (index_base, start_vol, stop_vol)
}

#[cfg(test)]
#[path = "mb_index_tests.rs"]
mod tests;
