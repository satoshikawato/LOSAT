//! blastdb_aliastool's `Run` (`CBlastDBAliasApp::Run`): the mode its arguments choose, the files
//! it writes, and its errors as `CATCH_ALL` reports them.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::{ArgException, Args};
use crate::ncbi_app::files::{self, Cin, LineEnds};
use crate::ncbi_app::{database_error, Diag};
use crate::ncbi_exception::NcbiException;
use crate::ncbi_str::string_to_int8;
use crate::seqdb::alias::{base_name, path_base_name};
use crate::seqdb::{id_list, seqidlist, SearchPath, SeqDb, SeqDbError, SeqType};
use crate::wording::{DATABASE_ERROR, QUERY_OPTIONS_ERROR, USER_MANUAL_FOOTER};
use crate::writedb::alias_file::{self, AliasError, FilterType, CLOCK_ERROR};
use crate::writedb::binary_list::gi_list_bytes;
use crate::writedb::seqidlist_writer::{seqidlist_bytes, SeqidlistError};
use crate::writedb::time::LocalTime;

use super::PROGRAM;

/// What ends `Run` before its end, as `CATCH_ALL` reports it (blast_app_util.hpp).
enum Failure {
    /// A `CSeqDBException` of the program or of the writer: `BLAST Database error: <text>`,
    /// exit status 2.
    SeqDb(Vec<u8>),
    /// The reader's error (`ncbi_app::database_error`).
    Reader(SeqDbError),
    /// A `CInputException`: `BLAST query/options error: <text>` and the manual's line, exit
    /// status 1.
    Input(&'static str),
    /// A `CArgException` thrown in `Run` (`AsInputFile`, an excluded value):
    /// `Command line argument error: <text>`, exit status 1.
    Arg(Vec<u8>),
    /// Another exception of NCBI's: `Error: <what()>`, exit status 255.
    Exception(NcbiException),
    /// LOSAT's refusal (a registry value LOSAT does not read) or error (the clock).
    Refused(String),
}

impl From<ArgException> for Failure {
    fn from(error: ArgException) -> Self {
        Failure::Arg(error.message().to_vec())
    }
}

impl From<AliasError> for Failure {
    fn from(error: AliasError) -> Self {
        match error {
            AliasError::SeqDb(text) => Failure::SeqDb(text.into_bytes()),
            AliasError::Reader(error) => Failure::Reader(error),
            AliasError::Clock => Failure::Refused(CLOCK_ERROR.to_string()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:628-665
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-255
/// `CBlastDBAliasApp::Run` on `args`; the exit status. Its files are written on the file system.
pub(super) fn run(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    match run_mode(args, settings, stdin, out, diag) {
        Ok(status) => status,
        Err(Failure::SeqDb(text)) => {
            diag.post_bytes(&[DATABASE_ERROR.as_bytes(), b" ", &text].concat());
            2
        }
        Err(Failure::Reader(error)) => database_error(diag, error),
        Err(Failure::Input(text)) => {
            diag.post(&format!("{QUERY_OPTIONS_ERROR} {text}"));
            diag.post(USER_MANUAL_FOOTER);
            1
        }
        Err(Failure::Arg(text)) => {
            diag.post_bytes(&[b"Command line argument error: ".as_slice(), &text].concat());
            1
        }
        Err(Failure::Exception(exception)) => {
            diag.post_bytes(&[b"Error: ".as_slice(), &exception.what()].concat());
            255
        }
        Err(Failure::Refused(reason)) => diag.refuse(&reason),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:103-115
/// The mode of `x_GetOperationMode` (its last test that holds wins), then its work.
fn run_mode(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    if args.given("seqid_file_info") {
        seqid_file_info(args, settings, out)?;
        Ok(0)
    } else if args.given("seqid_file_in") {
        convert_seqid_file(args, settings, stdin, diag)
    } else if args.given("gi_file_in") {
        convert_gi_file(args, stdin, diag)?;
        Ok(0)
    } else {
        create_alias_file(args, settings, stdin, diag)?;
        Ok(0)
    }
}

/// The database search path of NCBI's reader (`CSeqDBAtlas::GenerateSearchPath`: the current
/// directory, `BLASTDB`, the registry's `[BLAST] BLASTDB`), made where a mode first needs it.
fn search_path(settings: &ApplicationSettings) -> Result<SearchPath, Failure> {
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => {
            return Err(Failure::Refused(
                crate::algorithm::common::database::refusal(reason, PROGRAM).to_string(),
            ))
        }
    };
    SearchPath::from_environment(registry).map_err(Failure::Reader)
}

/// The bytes of the input file of the argument `name` (`AsInputFile`, `-` is `stdin`), as the
/// program's loop reads them; a file that does not open is the `CArgException` of `Run`.
fn input_file(args: &Args, name: &str, stdin: &mut Cin<'_>) -> Result<Vec<u8>, Failure> {
    let value = args.value(name).unwrap_or_default();
    let mut file = files::open_input_file(name, value, stdin)?;
    let mut bytes = Vec::new();
    // NCBI's loop ends at the stream's failure: the bytes read before it are the input.
    let _ = file.read_to_end(&mut bytes);
    Ok(bytes)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:345-384
// `NcbiGetlineEOL` reads the lines (on Unix a `\r` stays in the line, and `StringToInt8` of its
// default flags then fails); a line that does not convert is a warning, which blastdb_aliastool's
// post level does not show, and is skipped.
/// The GIs of a text GI list's bytes (`ConvertGiFile`): each line that is not empty and does not
/// start with `#`, as an `Int8`.
fn gi_ids(bytes: &[u8]) -> Vec<i64> {
    files::input_lines(bytes, LineEnds::HOST)
        .into_iter()
        .filter(|line| !line.starts_with(b"#"))
        .filter_map(|line| string_to_int8(line, 0).ok())
        .collect()
}

/// `CNcbiOfstream` of `path` in binary mode and the write of `bytes`; a failure is not reported,
/// as NCBI's stream does not report it.
fn write_file(path: &Path, bytes: &[u8]) {
    if let Ok(mut file) = std::fs::File::create(path) {
        let _ = file.write_all(bytes);
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:634-648
// `-gi_file_out` was made or emptied as the reader read it (`fPreOpen`); `Run` opens its name
// again. `CFile::Exists` is true for a regular file.
/// `-gi_file_in`: the binary GI list `-gi_file_out` (or `<gi_file_in>.bgl`) of the text list.
fn convert_gi_file(args: &Args, stdin: &mut Cin<'_>, diag: &mut Diag) -> Result<(), Failure> {
    let bytes = input_file(args, "gi_file_in", stdin)?;
    let gi_file_out: PathBuf = match args.os_value("gi_file_out") {
        Some(name) => PathBuf::from(name),
        None => PathBuf::from(format!(
            "{}.bgl",
            args.value("gi_file_in").unwrap_or_default()
        )),
    };
    let ids = gi_ids(&bytes);
    write_file(&gi_file_out, &gi_list_bytes(&ids));
    diag.post(&format!("Converted {} GIs into binary GI file", ids.len()));
    if !gi_file_out.is_file() {
        return Err(Failure::SeqDb(
            [
                gi_file_out.as_os_str().as_encoded_bytes(),
                b" not written".as_slice(),
            ]
            .concat(),
        ));
    }
    Ok(())
}

/// The sequence type of `-dbtype` (`prot`, else `nucl`), read as `Run` reads it: excluded by a
/// given argument, it is the reader's exception.
fn seq_type(args: &Args) -> Result<SeqType, Failure> {
    Ok(match args.as_string("dbtype")? {
        "prot" => SeqType::Protein,
        _ => SeqType::Nucleotide,
    })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3439-3470
// With the default flags an empty piece between two spaces, or before the first or after the
// last, is kept; an empty text has no piece.
/// `NStr::Split` of `text` at single spaces with its default flags.
fn split_names(text: &str) -> Vec<String> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(' ').map(str::to_string).collect()
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:548-570
// The file's lines are read by `getline` and their white space is cut (`NStr::TruncateSpaces`).
/// `x_GetDbsToAggregate`: the names of `dbs` split at spaces, else the non-empty lines of the
/// file `file`.
fn names_to_aggregate(
    args: &Args,
    dbs: &str,
    file: &str,
    stdin: &mut Cin<'_>,
) -> Result<Vec<String>, Failure> {
    if let Some(names) = args.value(dbs) {
        return Ok(split_names(names));
    }
    let bytes = input_file(args, file, stdin)?;
    Ok(bytes
        .split(|&byte| byte == b'\n')
        .map(|line| String::from_utf8_lossy(line.trim_ascii()).into_owned())
        .filter(|line| !line.is_empty())
        .collect())
}

/// Whether the argument `name` has a value (`HasValue`, also the `bool` of `args[name]`).
fn has(args: &Args, name: &str) -> bool {
    args.given(name)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:386-517
// The text GI list given to `-gilist` is converted first to `<out>.<p|n>.gil`, which the alias
// file names (`ConvertGiFile` with the names: its log line names both files). A `-dblist` is
// opened by the reader before the alias file is made (`CSeqDBExpert`), a `-dblist_file` is not.
// `-seqidlist`, `-taxidlist` and `-oid_masks` are ignored with `-dblist` or `-dblist_file`.
/// `CreateAliasFile`: the alias file that the arguments name, then the VDB names.
fn create_alias_file(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    diag: &mut Diag,
) -> Result<(), Failure> {
    let has_list = ["gilist", "seqidlist", "taxidlist", "oid_masks"]
        .iter()
        .any(|&name| has(args, name));
    if has(args, "db") && !has_list {
        return Err(Failure::Input(
            "Either gilist, seqid_list, taxidlist or oid_masks must be specified if database name \
             is used",
        ));
    }
    let title = match args.value("title") {
        Some(title) => title.to_string(),
        None => match args.value("db") {
            Some(db) => {
                let limit = if let Some(list) = args.value("gilist") {
                    list.to_string()
                } else if let Some(list) = args.value("seqidlist") {
                    list.to_string()
                } else if let Some(list) = args.value("taxidlist") {
                    list.to_string()
                } else if args.integer("oid_masks").is_some_and(exclude_model) {
                    "exclude model oid masks".to_string()
                } else {
                    String::new()
                };
                format!("{} limited by {limit}", path_base_name(db))
            }
            None => String::new(),
        },
    };
    let seqtype = seq_type(args)?;
    let out = args.value("out").unwrap_or_default();
    let mut gilist = args.value("gilist").unwrap_or_default().to_string();
    if !gilist.is_empty() {
        if !Path::new(&gilist).is_file() {
            return Err(Failure::SeqDb(format!("{gilist} not found").into_bytes()));
        }
        let bytes = std::fs::read(&gilist).map_err(|source| {
            Failure::Reader(SeqDbError::Open {
                path: PathBuf::from(&gilist),
                source,
            })
        })?;
        let binary =
            id_list::is_binary_gi_list(&bytes).map_err(|text| Failure::SeqDb(text.into_bytes()))?;
        if !binary {
            let mol_type = args.as_string("dbtype")?.chars().next().unwrap_or('p');
            let ifname = std::mem::replace(&mut gilist, format!("{out}.{mol_type}.gil"));
            let mut text = Vec::new();
            if let Ok(mut input) = files::open_text_file(&ifname) {
                let _ = input.read_to_end(&mut text);
            }
            let ids = gi_ids(&text);
            write_file(Path::new(&gilist), &gi_list_bytes(&ids));
            diag.post(&format!(
                "Converted {} GIs from {ifname} to binary format in {gilist}",
                ids.len()
            ));
        }
    }
    if let Some(dblist) = args.value("dblist") {
        SeqDb::open_in(dblist, seqtype, &search_path(settings)?).map_err(Failure::Reader)?;
    }
    let created = if has(args, "dblist") || has(args, "dblist_file") {
        let databases = names_to_aggregate(args, "dblist", "dblist_file", stdin)?;
        Some(alias_file::create_alias_file(
            out,
            &databases,
            seqtype,
            &gilist,
            &title,
            FilterType::GiList,
            &search_path(settings)?,
        )?)
    } else if let Some(num_volumes) = args.integer("num_volumes") {
        Some(alias_file::create_volumes_alias_file(
            out,
            u32::try_from(num_volumes).unwrap_or_default(),
            seqtype,
            &title,
            &search_path(settings)?,
        )?)
    } else if let Some(db) = args.value("db") {
        list_alias_file(args, db, out, seqtype, &gilist, &title, settings)?
    } else {
        None
    };
    if let Some(line) = created {
        diag.post(&line);
    }
    if has(args, "vdblist") || has(args, "vdblist_file") {
        let (extension, append) = match (seqtype, has(args, "dblist") || has(args, "dblist_file")) {
            (SeqType::Protein, true) => (".pal", true),
            (SeqType::Nucleotide, true) => (".nal", true),
            (SeqType::Protein, false) => (".pvl", false),
            (SeqType::Nucleotide, false) => (".nvl", false),
        };
        let title = match append {
            true => "",
            false => args.value("title").unwrap_or_default(),
        };
        add_vdbs(args, &format!("{out}{extension}"), append, title, stdin)?;
    }
    Ok(())
}

/// Whether `-oid_masks` names the exclude-model masks (`fExcludeModel` of NCBI's `int`).
fn exclude_model(mask: i32) -> bool {
    (mask as u64) & crate::seqdb::oid_mask::EXCLUDE_MODEL != 0
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:461-494
// `-db` is split at spaces as `-dblist` is (the overload of `CWriteDB_CreateAliasFile` with a
// name, writedb.cpp:330-340). A `-seqidlist` or `-taxidlist` file must exist; the OID mask
// database writes no log line.
/// The alias file of `-db` and its list or OID masks; `Ok(Some)` is the log line.
fn list_alias_file(
    args: &Args,
    db: &str,
    out: &str,
    seqtype: SeqType,
    gilist: &str,
    title: &str,
    settings: &ApplicationSettings,
) -> Result<Option<String>, Failure> {
    let (list, filter) = if has(args, "gilist") {
        (gilist, FilterType::GiList)
    } else if let Some(list) = args.value("seqidlist") {
        (list, FilterType::SeqIdList)
    } else if let Some(list) = args.value("taxidlist") {
        (list, FilterType::TaxIdList)
    } else if let Some(mask) = args.integer("oid_masks") {
        alias_file::create_oid_mask_db(db, out, seqtype, mask, title, &search_path(settings)?)?;
        return Ok(None);
    } else {
        return Ok(None);
    };
    if filter != FilterType::GiList && !list.is_empty() && !Path::new(list).is_file() {
        return Err(Failure::SeqDb(format!("{list} not found").into_bytes()));
    }
    Ok(Some(alias_file::create_alias_file(
        out,
        &split_names(db),
        seqtype,
        list,
        title,
        filter,
        &search_path(settings)?,
    )?))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:519-546
// An empty list posts a warning, which blastdb_aliastool's post level does not show, and writes
// nothing. A file that does not open drops the writes.
/// `x_AddVDBsToAliasFile`: the VDB names appended to the alias file `filename` (`append`), or a
/// new VDB alias file with its header and `title` (none when empty).
fn add_vdbs(
    args: &Args,
    filename: &str,
    append: bool,
    title: &str,
    stdin: &mut Cin<'_>,
) -> Result<(), Failure> {
    let vdbs = names_to_aggregate(args, "vdblist", "vdblist_file", stdin)?;
    if vdbs.is_empty() {
        return Ok(());
    }
    let mut text = String::new();
    if !append {
        let now = LocalTime::now().ok_or_else(|| Failure::Refused(CLOCK_ERROR.to_string()))?;
        text.push_str(&format!("#\n# Alias file created {}\n#\n", now.log_time()));
    }
    if !title.is_empty() {
        text.push_str(&format!("TITLE {title}\n"));
    }
    text.push_str("VDBLIST ");
    for name in &vdbs {
        text.push_str(&format!("\"{name}\" "));
    }
    text.push('\n');
    let file = std::fs::OpenOptions::new()
        .write(true)
        .create(true)
        .append(append)
        .truncate(!append)
        .open(filename);
    if let Ok(mut file) = file {
        let _ = file.write_all(text.as_bytes());
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:572-618
// The output file is made before the input is read and checked: a list that fails leaves it
// empty. The default title is the input's name after its last path separator
// (`CSeqDB_Path::FindFileName`). With `-seqid_db` and no `-seqid_dbtype` the reader tries a
// protein database, then a nucleotide one (`CSeqDB` of `eUnknown`).
/// `-seqid_file_in`: the binary seqid list `-seqid_file_out` (or `<seqid_file_in>.bsl`); the exit
/// status (255 for a list without ids, which NCBI's writer returns as -1).
fn convert_seqid_file(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    diag: &mut Diag,
) -> Result<i32, Failure> {
    let input_name = args.value("seqid_file_in").unwrap_or_default();
    let mut input = files::open_input_file("seqid_file_in", input_name, stdin)?;
    let out_filename = match args.value("seqid_file_out") {
        Some(name) => name.to_string(),
        None => format!("{input_name}.bsl"),
    };
    let title = match args.value("seqid_title") {
        Some(title) => title.to_string(),
        None => base_name(input_name).to_string(),
    };
    let mut output = std::fs::File::create(&out_filename).ok();
    // The input is read after the output is made: an output of the input's name leaves nothing
    // to read.
    let mut bytes = Vec::new();
    let _ = input.read_to_end(&mut bytes);
    let seqid_list: Vec<Vec<u8>> = files::input_lines(&bytes, LineEnds::HOST)
        .into_iter()
        .filter(|line| !line.starts_with(b"#"))
        .map(<[u8]>::to_vec)
        .collect();
    let db = match args.value("seqid_db") {
        Some(name) => {
            let path = search_path(settings)?;
            let types: &[SeqType] = match args.value("seqid_dbtype") {
                Some(dbtype) if dbtype.starts_with('p') => &[SeqType::Protein],
                Some(_) => &[SeqType::Nucleotide],
                None => &[SeqType::Protein, SeqType::Nucleotide],
            };
            Some(open_guessed(name, types, &path).map_err(Failure::Reader)?)
        }
        None => None,
    };
    match seqidlist_bytes(&seqid_list, title.as_bytes(), db.as_ref()) {
        Ok(list) => {
            if let Some(output) = &mut output {
                let _ = output.write_all(&list);
            }
            Ok(0)
        }
        Err(SeqidlistError::SeqDb(text)) => Err(Failure::SeqDb(text)),
        Err(SeqidlistError::InvalidId { posted, message }) => {
            diag.err_post_bytes(&posted);
            Err(Failure::SeqDb(message))
        }
        Err(SeqidlistError::Empty) => {
            diag.err_post("Empty seqid list");
            Ok(255)
        }
        Err(SeqidlistError::Reader(error)) => Err(Failure::Reader(error)),
        Err(SeqidlistError::Clock) => Err(Failure::Refused(CLOCK_ERROR.to_string())),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:110-153
// A protein database that does not open as one (a `CSeqDBException`) is tried as a nucleotide
// one; that second error is reported.
/// `CSeqDB` of `name` and the types `types` in order (both for `eUnknown`).
fn open_guessed(name: &str, types: &[SeqType], path: &SearchPath) -> Result<SeqDb, SeqDbError> {
    let mut types = types.iter().peekable();
    loop {
        let seqtype = *types.next().expect("at least one sequence type is tried");
        match SeqDb::open_in(name, seqtype, path) {
            Err(error) if types.peek().is_some() && error.is_seqdb_exception() => {}
            opened => return opened,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:620-626
/// `-seqid_file_info`: what the binary seqid list on the search path says of itself, on `out`.
fn seqid_file_info(
    args: &Args,
    settings: &ApplicationSettings,
    out: &mut dyn Write,
) -> Result<(), Failure> {
    let name = args.value("seqid_file_info").unwrap_or_default();
    let bytes =
        seqidlist::read_resolved(name, &search_path(settings)?).map_err(Failure::Exception)?;
    let text = seqidlist::info_text(&bytes).map_err(|text| Failure::SeqDb(text.into_bytes()))?;
    let _ = out.write_all(&text);
    let _ = out.flush();
    Ok(())
}

#[cfg(test)]
#[path = "work_tests.rs"]
mod tests;
