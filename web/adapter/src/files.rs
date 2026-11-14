//! Named input files (docs/web/abi_v2.md §4 `file_begin`): files that an argv names by a value
//! of an option, which the host gives in chunks under that name. Kind 1 is a WindowMasker unit
//! counts file: the value of BLASTN's `-window_masker_db`, or the file of `-window_masker_taxid
//! <taxid>`, `windowmasker/<taxid>/wmasker.obinary` (`WINDOW_MASKER_PATH`). Kind 2 is a registry
//! file of NCBI's application start: the file of `-conffile`, which the start of a run reads
//! (`registry_files`; line RG item E1w-3a).
//!
//! The command line opens the counts file by its name twice: NCBI's check of the format
//! (`DiscoverStatType`) and its reader (`CSeqMaskerIstatFactory::create`). A registered file
//! answers both from what was read once at `file_end`: the bytes are checked and parsed there,
//! kept as the parsed counts, and dropped. A run of the argv gives the engine a source of the
//! registered files (`source`), which answers a name that is not registered as the command line
//! answers a file that does not exist.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::Cursor;
use std::sync::{Arc, Mutex, OnceLock};

use losat::web_adapter::{
    discover_stat_type_in, CountsSource, MaskerParams, StatFile, StatType, UnitCounts, WinMaskError,
};
use sha2::{Digest, Sha256};

use crate::json;

/// The kind of a WindowMasker unit counts file (`-window_masker_db`, `-window_masker_taxid`).
pub const KIND_WINDOW_MASKER_COUNTS: u32 = 1;

/// The kind of a registry file of NCBI's application start (`-conffile`).
pub const KIND_REGISTRY: u32 = 2;

/// The directory of BLASTN's `-window_masker_taxid` in LOSAT Web, where the command line reads
/// NCBI's registry entry `[WINDOW_MASKER] WINDOW_MASKER_PATH`: the counts file of a taxid is the
/// file registered as `windowmasker/<taxid>/wmasker.obinary` (the engine joins the names as NCBI
/// does, `taxid_stat_file`), and a taxid without one is the command line's taxid without a file.
pub const WINDOW_MASKER_PATH: &str = "windowmasker";

/// A file whose chunks are coming.
struct Loading {
    kind: u32,
    name: String,
    size: usize,
    bytes: Vec<u8>,
    sha: Sha256,
}

/// A registered counts file: what NCBI's check of the format found, and what its reader read, as
/// `Result`s, so that a run reports them where the command line does. The reader's result is kept
/// for every format: `-window_masker_db` reads the file only after the check accepted one of the
/// two binary formats (`blast_args.cpp:446-458`), but the file of `-window_masker_taxid` is read
/// without that check (`WindowMaskerTaxidToDb`), and the reader's own errors then come with the
/// search ("unrecognized unit counts format" for an unknown format).
struct Counts {
    name: String,
    sha256: String,
    found: Result<StatFile, WinMaskError>,
    counts: Result<Arc<UnitCounts>, WinMaskError>,
}

/// A registered registry file: its bytes, which the start of each run reads (`registry_files`).
struct RegistryFile {
    name: String,
    bytes: Vec<u8>,
}

#[derive(Default)]
struct Files {
    next: u32,
    loading: HashMap<u32, Loading>,
    registered: HashMap<u32, Arc<Counts>>,
    registry: HashMap<u32, Arc<RegistryFile>>,
}

impl Files {
    fn names(&self) -> impl Iterator<Item = &str> {
        self.loading
            .values()
            .map(|file| file.name.as_str())
            .chain(self.registered.values().map(|file| file.name.as_str()))
            .chain(self.registry.values().map(|file| file.name.as_str()))
    }
}

fn files() -> &'static Mutex<Files> {
    static FILES: OnceLock<Mutex<Files>> = OnceLock::new();
    FILES.get_or_init(Mutex::default)
}

/// `file_begin`: a handle for the chunks of a file of `kind` named `name`, of `size` bytes. The
/// memory for the bytes is taken here, so that a file that does not fit fails before its chunks
/// are copied.
pub fn begin(kind: u32, name: &str, size: usize) -> Result<u32, String> {
    let what = match kind {
        KIND_WINDOW_MASKER_COUNTS => "the WindowMasker counts",
        KIND_REGISTRY => "the registry file",
        _ => return Err(format!("unknown input file kind {kind}")),
    };
    if name.is_empty() {
        return Err("an input file needs a name".to_string());
    }
    let mut files = files().lock().expect("files");
    if files.names().any(|known| known == name) {
        return Err(format!(
            "an input file named {name} is registered; release it before registering another"
        ));
    }
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(size)
        .map_err(|_| format!("not enough memory for {what} {name} ({size} bytes)"))?;
    let handle = files
        .next
        .checked_add(1)
        .filter(|handle| i32::try_from(*handle).is_ok())
        .ok_or_else(|| "input file handle space exhausted".to_string())?;
    files.next = handle;
    files.loading.insert(
        handle,
        Loading {
            kind,
            name: name.to_string(),
            size,
            bytes,
            sha: Sha256::new(),
        },
    );
    Ok(handle)
}

/// `file_chunk`: the next bytes of the file, in order.
pub fn chunk(file: u32, data: &[u8]) -> Result<(), String> {
    let mut files = files().lock().expect("files");
    let loading = files
        .loading
        .get_mut(&file)
        .ok_or_else(|| format!("unknown input file {file}"))?;
    if data.len() > loading.size - loading.bytes.len() {
        let message = format!(
            "the input file {} has more than its {} bytes",
            loading.name, loading.size
        );
        files.loading.remove(&file);
        return Err(message);
    }
    loading.bytes.extend_from_slice(data);
    loading.sha.update(data);
    Ok(())
}

/// `file_end`: checks and reads the file as the command line would, keeps the results, drops the
/// bytes, and returns the *file* JSON. An error of the file's contents does not fail here: the
/// run reports it as the command line does.
pub fn end(file: u32) -> Result<String, String> {
    let loading = files()
        .lock()
        .expect("files")
        .loading
        .remove(&file)
        .ok_or_else(|| format!("unknown input file {file}"))?;
    if loading.bytes.len() != loading.size {
        return Err(format!(
            "the input file {} ended after {} of its {} bytes",
            loading.name,
            loading.bytes.len(),
            loading.size
        ));
    }
    let Loading {
        kind,
        name,
        size,
        bytes,
        sha,
    } = loading;
    let mut sha256 = String::with_capacity(64);
    for byte in sha.finalize() {
        let _ = write!(sha256, "{byte:02x}");
    }
    if kind == KIND_REGISTRY {
        // The bytes are kept: the start of each run reads the file as the command line's start
        // reads it from the disk (`registry_files`), its errors included.
        let mut out = String::new();
        let _ = write!(out, "{{\"handle\":{file},\"name\":");
        json::string(&mut out, &name);
        let _ = write!(
            out,
            ",\"kind\":\"registry\",\"size\":{size},\"sha256\":\"{sha256}\"}}"
        );
        files()
            .lock()
            .expect("files")
            .registry
            .insert(file, Arc::new(RegistryFile { name, bytes }));
        return Ok(out);
    }
    let found = discover_stat_type_in(&mut Cursor::new(bytes.as_slice()));
    // The command line's reader with blastn's settings (`s_BuildSeqMasker`); for a format other
    // than the two binary ones it reads the format alone and returns its error.
    let counts = UnitCounts::read(&mut Cursor::new(bytes.as_slice()), &MaskerParams::default())
        .map(Arc::new);
    drop(bytes);
    let counts = Arc::new(Counts {
        name,
        sha256,
        found,
        counts,
    });
    let response = file_json(file, size, &counts);
    files()
        .lock()
        .expect("files")
        .registered
        .insert(file, counts);
    Ok(response)
}

/// `file_release`: forgets the file (registered or still loading) and frees its memory.
pub fn release(file: u32) -> Result<(), String> {
    let mut files = files().lock().expect("files");
    let removed = files.registered.remove(&file).is_some()
        || files.registry.remove(&file).is_some()
        || files.loading.remove(&file).is_some();
    if removed {
        Ok(())
    } else {
        Err(format!("unknown input file {file}"))
    }
}

/// The two formats that NCBI's check of `-window_masker_db` accepts.
fn is_binary(stat_type: StatType) -> bool {
    matches!(stat_type, StatType::OBinary | StatType::Binary)
}

fn stat_type_name(stat_type: StatType) -> &'static str {
    match stat_type {
        StatType::Ascii => "ascii",
        StatType::Binary => "binary",
        StatType::OAscii => "oascii",
        StatType::OBinary => "obinary",
        StatType::Unknown => "unknown",
    }
}

/// The *file* JSON (docs/web/abi_v2.md §9).
fn file_json(handle: u32, size: usize, file: &Counts) -> String {
    let mut out = String::new();
    let _ = write!(out, "{{\"handle\":{handle},\"name\":");
    json::string(&mut out, &file.name);
    let _ = write!(
        out,
        ",\"kind\":\"windowmasker_counts\",\"size\":{size},\"sha256\":\"{}\",\"stat_type\":",
        file.sha256
    );
    json::optional_string(
        &mut out,
        file.found
            .as_ref()
            .ok()
            .map(|stat| stat_type_name(stat.stat_type)),
    );
    out.push_str(",\"unit_size\":");
    json::optional(
        &mut out,
        file.counts.as_ref().ok().map(|counts| counts.unit_size()),
    );
    out.push_str(",\"error\":");
    // The reader's error of a file that the check refused is the check's: only the check's is told.
    let error = match (&file.found, &file.counts) {
        (Err(error), _) => Some(error.to_string()),
        (Ok(stat), Err(error)) if is_binary(stat.stat_type) => Some(error.to_string()),
        _ => None,
    };
    json::optional_string(&mut out, error.as_deref());
    out.push('}');
    out
}

/// The registered counts files of a run, by name. The engine reads them through this source
/// instead of the file system; the table is copied when the run starts, so the search's threads
/// never take the lock of the registrations.
#[derive(Default)]
pub struct RegisteredCounts {
    files: HashMap<String, Arc<Counts>>,
}

impl std::fmt::Debug for RegisteredCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut list = f.debug_map();
        for (name, file) in &self.files {
            list.entry(name, &file.sha256);
        }
        list.finish()
    }
}

impl RegisteredCounts {
    fn file(&self, name: &str) -> Result<&Counts, WinMaskError> {
        self.files
            .get(name)
            .map(|file| &**file)
            .ok_or_else(|| WinMaskError::Open {
                name: name.to_string(),
            })
    }
}

impl CountsSource for RegisteredCounts {
    fn discover(&self, name: &str) -> Result<StatFile, WinMaskError> {
        self.file(name)?.found.clone()
    }

    fn counts(&self, name: &str, params: &MaskerParams) -> Result<Arc<UnitCounts>, WinMaskError> {
        let file = self.file(name)?;
        if *params != MaskerParams::default() {
            return Err(WinMaskError::CreateFail {
                message: format!(
                    "LOSAT Web holds the WindowMasker counts {name} read with BLASTN's settings only"
                ),
            });
        }
        file.counts.clone()
    }

    /// A registered file is a file of its name (compared byte for byte, as the names of
    /// `-window_masker_db`); nothing else is.
    fn is_file(&self, path: &std::path::Path) -> bool {
        path.to_str()
            .is_some_and(|name| self.files.contains_key(name))
    }
}

/// Whether a counts file of this name is registered (`file_end` succeeded and it was not
/// released).
pub fn registered(name: &str) -> bool {
    files()
        .lock()
        .expect("files")
        .registered
        .values()
        .any(|file| file.name == name)
}

/// The source of the counts files for a run: the files registered now, by name.
pub fn source() -> Arc<dyn CountsSource> {
    let files = files().lock().expect("files");
    Arc::new(RegisteredCounts {
        files: files
            .registered
            .values()
            .map(|file| (file.name.clone(), file.clone()))
            .collect(),
    })
}

/// The registry files of a run (kind 2), by name: the files that its application start reads
/// instead of the file system (`losat::web_adapter::start_run`). A registered file is the file of
/// its name in the command line's working directory, so a leading `./` (the working directory on
/// NCBI's registry search path, or written in the name) names the same file (`working_name`);
/// nothing else is a file. The table is copied when the run starts.
pub struct RegisteredRegistry {
    files: HashMap<String, Arc<RegistryFile>>,
}

/// The name of `path` in the working directory: without its leading `.` components (`./c.ini` and
/// `c.ini` are one file); `None` for a path that is not UTF-8 or names the directory itself.
fn working_name(path: &std::path::Path) -> Option<String> {
    let rest: std::path::PathBuf = path
        .components()
        .skip_while(|part| matches!(part, std::path::Component::CurDir))
        .collect();
    let name = rest.to_str()?;
    (!name.is_empty()).then(|| name.to_string())
}

impl RegisteredRegistry {
    fn file(&self, path: &std::path::Path) -> Option<&RegistryFile> {
        working_name(path)
            .and_then(|name| self.files.get(&name))
            .map(|file| &**file)
    }
}

impl losat::web_adapter::RegistryFiles for RegisteredRegistry {
    fn is_file(&self, path: &std::path::Path) -> bool {
        self.file(path).is_some()
    }

    fn opens(&self, path: &std::path::Path) -> bool {
        self.file(path).is_some()
    }

    fn read(&self, path: &std::path::Path) -> Option<std::io::Result<Vec<u8>>> {
        self.file(path).map(|file| Ok(file.bytes.clone()))
    }
}

/// The registry files for a run: the files of kind 2 registered now.
pub fn registry_files() -> RegisteredRegistry {
    let files = files().lock().expect("files");
    RegisteredRegistry {
        files: files
            .registry
            .values()
            .filter_map(|file| {
                working_name(std::path::Path::new(&file.name)).map(|name| (name, file.clone()))
            })
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts_bytes(name: &str) -> Vec<u8> {
        std::fs::read(format!(
            "{}/../../LOSAT/tests/winmask/counts/{name}",
            env!("CARGO_MANIFEST_DIR")
        ))
        .expect("a counts fixture")
    }

    fn sha256_hex(bytes: &[u8]) -> String {
        Sha256::digest(bytes)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }

    fn register(name: &str, bytes: &[u8], chunk_size: usize) -> (u32, String) {
        let file = begin(KIND_WINDOW_MASKER_COUNTS, name, bytes.len()).expect("begin");
        for part in bytes.chunks(chunk_size.max(1)) {
            chunk(file, part).expect("chunk");
        }
        (file, end(file).expect("end"))
    }

    // The tests share the process's one table, so each uses names of its own.

    #[test]
    fn chunk_sizes_give_the_same_hash_and_counts() {
        let bytes = counts_bytes("mg1655_u9.obinary");
        let whole = UnitCounts::read(&mut Cursor::new(bytes.as_slice()), &MaskerParams::default())
            .expect("the fixture reads");
        for (index, size) in [1usize << 20, 4096, 7, bytes.len()].into_iter().enumerate() {
            let name = format!("chunks{index}.obinary");
            let (file, response) = register(&name, &bytes, size);
            assert!(response.contains(&format!("\"sha256\":\"{}\"", sha256_hex(&bytes))));
            assert!(response.contains("\"stat_type\":\"obinary\",\"unit_size\":9,\"error\":null"));
            let source = source();
            let counts = source
                .counts(&name, &MaskerParams::default())
                .expect("registered counts");
            assert_eq!(counts.unit_size(), whole.unit_size());
            assert_eq!(counts.threshold(), whole.threshold());
            for unit in [0u32, 1, 12345, 0x3ffff] {
                assert_eq!(counts.true_count(unit), whole.true_count(unit));
            }
            release(file).expect("release");
        }
    }

    // The file of `-window_masker_taxid` is read without NCBI's check of the format: the counts of a
    // registered file are the command line's reader's result for any format.
    #[test]
    fn the_counts_of_any_format_are_the_readers() {
        for (name, bytes) in [
            ("empty.any", &b""[..]),
            ("short.any", &b"abc"[..]),
            ("text.any", &b"9\nAAAAAAAAA 3\n"[..]),
        ] {
            let (file, response) = register(name, bytes, 4);
            let want = UnitCounts::read(&mut Cursor::new(bytes), &MaskerParams::default())
                .map(|_| ())
                .expect_err("not counts LOSAT reads");
            let got = source()
                .counts(name, &MaskerParams::default())
                .map(|_| ())
                .expect_err("the reader's error");
            assert_eq!(got, want, "{name}");
            assert!(response.contains("\"unit_size\":null"), "{response}");
            release(file).expect("release");
        }
    }

    #[test]
    fn registered_names_are_the_files_of_the_taxid_directory() {
        use std::path::Path;
        let bytes = counts_bytes("mg1655_u9.obinary");
        let name = format!("{WINDOW_MASKER_PATH}/77/wmasker.obinary");
        let (file, _) = register(&name, &bytes, 1 << 16);
        let held = source();
        assert!(held.is_file(Path::new(&name)));
        for other in [
            "windowmasker/78/wmasker.obinary",
            "windowmasker//77/wmasker.obinary",
            "windowmasker/77",
        ] {
            assert!(!held.is_file(Path::new(other)), "{other}");
        }
        release(file).expect("release");
        assert!(held.is_file(Path::new(&name)), "a run keeps its table");
        assert!(!source().is_file(Path::new(&name)), "released");
    }

    #[test]
    fn sizes_are_kept() {
        let file = begin(KIND_WINDOW_MASKER_COUNTS, "size.over", 4).expect("begin");
        chunk(file, b"abc").expect("within the size");
        let error = chunk(file, b"de").expect_err("past the size");
        assert_eq!(error, "the input file size.over has more than its 4 bytes");
        assert!(end(file).is_err(), "the file was dropped");
        let file = begin(KIND_WINDOW_MASKER_COUNTS, "size.over", 4).expect("the name is free");
        chunk(file, b"abc").expect("chunk");
        assert_eq!(
            end(file).expect_err("short"),
            "the input file size.over ended after 3 of its 4 bytes"
        );
        assert!(release(file).is_err(), "a failed end registers nothing");
    }

    #[test]
    fn names_are_registered_once() {
        let bytes = counts_bytes("dm2L_u9.binary");
        let (file, response) = register("twice.binary", &bytes, 1 << 16);
        assert!(response.contains("\"stat_type\":\"binary\""));
        assert!(begin(KIND_WINDOW_MASKER_COUNTS, "twice.binary", 1).is_err());
        // One name is one file, of one kind.
        assert!(begin(KIND_REGISTRY, "twice.binary", 1).is_err());
        assert_eq!(
            begin(3, "other", 1).unwrap_err(),
            "unknown input file kind 3"
        );
        assert!(begin(KIND_WINDOW_MASKER_COUNTS, "", 1).is_err());
        release(file).expect("release");
        assert!(release(file).is_err());
        let (again, _) = register("twice.binary", &bytes, 1 << 16);
        assert_ne!(again, file);
        release(again).expect("release");
    }

    #[test]
    fn content_errors_register_and_come_with_the_run() {
        let bytes = counts_bytes("mg1655_u9.obinary");
        let (file, response) = register("short.obinary", &bytes[..bytes.len() / 2], 1 << 16);
        assert!(
            response.contains(
                "\"stat_type\":\"obinary\",\"unit_size\":null,\"error\":\"not enough data"
            ),
            "{response}"
        );
        let source = source();
        assert!(source.discover("short.obinary").is_ok());
        assert!(matches!(
            source.counts("short.obinary", &MaskerParams::default()),
            Err(WinMaskError::Format { .. })
        ));
        release(file).expect("release");

        let (file, response) = register("text.counts", b"9\nAAAAAAAAA 3\n", 4);
        assert!(response.contains("\"error\":null"), "{response}");
        assert!(
            !response.contains("\"stat_type\":\"obinary\""),
            "{response}"
        );
        release(file).expect("release");

        assert_eq!(
            source.discover("never.registered"),
            Err(WinMaskError::Open {
                name: "never.registered".to_string()
            })
        );
    }
}
