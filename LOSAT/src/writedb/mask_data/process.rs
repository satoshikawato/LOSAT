//! `x_ProcessMaskData` of makeblastdb (`c++/src/app/blastdb/makeblastdb.cpp`): the options
//! `-mask_data`, `-mask_id` and `-mask_desc` are split into lists, every mask
//! file is read with [`MaskFileReader`], and the masks are registered with the `BuildDatabase`
//! and collected for the sequences it will write.
//!
//! # The order in which `losat makeblastdb` calls things
//!
//! `makeblastdb::build_database_with_log` calls [`process_mask_data`] where `makeblastdb`'s
//! `x_BuildDatabase` does, in this order
//! (makeblastdb.cpp:1231-1296; `Run`, makeblastdb.cpp:1322-1331):
//!
//! ```text
//! 1. BuildDatabase::new(..)          CBuildDatabase constructor: the head of the log
//!                                    ("Building a new DB, current time", "New DB name:", ...)
//!                                    to the log stream (standard output, or -logfile)
//! 2. set_max_file_size(bytes)        after "Maximum file size: " (log, standard output)
//! 3. process_mask_data(..)           x_ProcessMaskData: one "Mask file: <name>" log line per file
//!                                    read; "Ignoring mask file ..." to the standard error stream
//! 4. add_fasta(..)                   x_ProcessInputData: the checks of the input file ("File X
//!                                    does not exist", "is empty") come after the mask files
//! 5. end_build(..)                   m_DB->EndBuild(): the closing lines of the log
//! 6. (CATCH_ALL)                     a failure of 1-5: the message and the status of
//!                                    `MaskDataError::into_makedb_error` /
//!                                    `MakedbError::exit_code`, one line on the standard error
//! 7. post_destruction_diagnostics    ~CBuildDatabase, when the application object goes: after
//!                                    the exception message and after the metadata file
//! ```
//!
//! A failure of step 3 ends the run at once: steps 4 and 5 are not run (no database is written
//! and the input file is not looked at), step 6 prints the message, and step 7 still posts the
//! "No sequences matched any of the masks provided." message when a mask list with masks was
//! read before the failure (`BuildDatabase::mask_data_source_mut` is set). A file that ends where
//! `more TRUE` announced a continuation ([`MaskFileError::EofInContinuation`]) is such a failure
//! with no message and exit status 0: makeblastdb leaves without a database and without an error
//! (`CATCH_ALL` ignores a `CIOException` that is not a flush failure).
//!
//! Streams: the log (steps 1, 2, 3 "Mask file:", 5) is makeblastdb's standard output; the
//! diagnostics (`ERR_POST`: "Ignoring mask file", the warnings of the FASTA reader, the message of
//! step 7) and the message of step 6 (`LOG_POST`) are on the standard error stream, in the order
//! of the events.
//!
//! The options `-mask_id` and `-mask_desc` depend on others in the argument parser
//! (makeblastdb.cpp:276-288): `-mask_id` needs `-mask_data`, `-mask_desc` needs `-mask_id`.
//! `-gi_mask_name` (which needs `-gi_mask`) replaces the file name in the registration
//! (makeblastdb.cpp:960-971); LOSAT writes no GI masks yet and refuses `-gi_mask` and
//! `-gi_mask_name` in the CLI, so that branch is not here.

use std::fmt;
use std::path::{Path, PathBuf};

use super::MaskedRangeSet;
use crate::ncbi_str::{string_to_int, truncate_spaces_str};
use crate::writedb::build_db::BuildDatabase;
use crate::writedb::fasta_source::Warn;
use crate::writedb::makeblastdb::MakedbError;
use crate::writedb::mask_info::{MaskFileError, MaskFileReader, MaskPost};
use crate::writedb::WriteDbError;

/// The values of the options that `x_ProcessMaskData` reads: `None` for an option that was not
/// given (`CArgValue::HasValue()` false), `Some("")` for one given with an empty value.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct MaskDataOptions<'a> {
    /// `-mask_data`.
    pub(crate) mask_data: Option<&'a str>,
    /// `-mask_id`.
    pub(crate) mask_id: Option<&'a str>,
    /// `-mask_desc`.
    pub(crate) mask_desc: Option<&'a str>,
    /// The directory in which a relative mask file name is looked up (`None`: the current
    /// directory, as makeblastdb does). The names in the log, the messages and the registry stay
    /// as given. A hook for tests, which cannot change the directory of the process.
    pub(crate) root: Option<&'a Path>,
}

/// Why `x_ProcessMaskData` stopped. Every one of them ends the run of makeblastdb (the exception
/// leaves `x_BuildDatabase`); [`MaskDataError::into_makedb_error`] gives the message and the exit
/// status.
#[derive(Debug)]
pub(crate) enum MaskDataError {
    /// An error of `x_ProcessMaskData` itself or of the database: `CInvalidDataException`
    /// ([`WriteDbError::Options`]), `CWriteDBException` ([`WriteDbError::Creation`]) of the
    /// registry, LOSAT's refusal of a negative algorithm id ([`WriteDbError::Unsupported`]).
    Ncbi(WriteDbError),
    /// A failure of the mask file reader: [`MaskFileError::Ncbi`] (the format guess and what
    /// NCBI's serial reader throws, NCBI's text), [`MaskFileError::Unsupported`],
    /// [`MaskFileError::EofInContinuation`].
    File(MaskFileError),
}

impl From<MaskFileError> for MaskDataError {
    fn from(error: MaskFileError) -> Self {
        MaskDataError::File(error)
    }
}

impl From<WriteDbError> for MaskDataError {
    fn from(error: WriteDbError) -> Self {
        MaskDataError::Ncbi(error)
    }
}

impl fmt::Display for MaskDataError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MaskDataError::Ncbi(error) => fmt::Display::fmt(error, f),
            MaskDataError::File(error) => fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for MaskDataError {}

impl MaskDataError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-171
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:212-215
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:236-247
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:256-263
    /// Where each outcome of `x_ProcessMaskData` becomes makeblastdb's (exit status, message);
    /// the only place that decides it. `None` is a run that ends with status 0 and no message.
    ///
    /// * [`MaskDataError::Ncbi`], [`MaskFileError::Ncbi`]: `CATCH_ALL` as NCBI (`BLAST options
    ///   error: ...`, status 1; `BLAST Database creation error: ...`, status 1; the message of a
    ///   failure of the file system, "Error: ...", status 255); what NCBI's serial reader throws
    ///   for a file it refuses ("Error: NCBI C++ Exception: ..." of a `CException`, status 255;
    ///   "BLAST ran out of memory" of a `std::bad_alloc`, status 4; "Error: ..." of the
    ///   `std::length_error` of a string longer than `std::string` holds, status 255); LOSAT's
    ///   refusal of a negative algorithm id ("losat makeblastdb: ...", status 1).
    /// * [`MaskFileError::Unsupported`]: NCBI reads the file and LOSAT does not
    ///   (`MakedbError::Unsupported`, status 1).
    /// * [`MaskFileError::EofInContinuation`]: `None`, exactly as NCBI (status 0, nothing is
    ///   printed here; the destructor's message follows).
    pub(crate) fn into_makedb_error(self) -> Option<MakedbError> {
        match self {
            MaskDataError::Ncbi(error) | MaskDataError::File(MaskFileError::Ncbi(error)) => {
                Some(MakedbError::Build(error))
            }
            MaskDataError::File(MaskFileError::EofInContinuation { .. }) => None,
            MaskDataError::File(MaskFileError::Unsupported(what)) => Some(
                MakedbError::Unsupported(format!("mask data that LOSAT does not read: {what}")),
            ),
        }
    }
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-275,282-301
/// `NStr::Split(text, ",", list)` (no flags): the empty string gives no token, any other string
/// one token for each comma plus one, empty tokens included ("," gives two, "a," gives "a" and
/// ""). The tokens are not trimmed.
fn split_commas(text: &str) -> Vec<&str> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(',').collect()
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:897-909
// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:3611-3616,3663-3668
/// `NStr::StringToNumeric<Int4>(text, &nid, fConvErr_NoThrow, 10)`: whether `text` converts
/// (`StringToInt` with no flag, `ncbi_str::string_to_int`; with `fConvErr_NoThrow` a failed
/// conversion gives 0 and sets errno).
fn is_numeric_int4(text: &str) -> bool {
    string_to_int(text.as_bytes(), 0).is_ok()
}

fn options_error(message: &str) -> MaskDataError {
    MaskDataError::Ncbi(WriteDbError::Options(message.to_string()))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:876-1000
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp (c++/include/corelib/ncbifile.hpp:4038-4042)
// `CFile::Exists` is true for a regular file (a link is followed); a directory, a device and a
// FIFO "do not exist" (the oracle prints "Ignoring mask file '/dev/null' as it does not
// exist."). A file that exists and cannot be opened is a failed stream, which the format guess
// finds neither binary nor text: "Unknown encoding for mask data in '<file>'" (the oracle, a
// file with mode 000).
/// `x_ProcessMaskData`. `db` is the `CBuildDatabase` (`m_DB`; the log `m_LogFile` is its log) and
/// `warn` takes the lines that NCBI posts with `ERR_POST` (without the final newline).
///
/// Returns at once, without a trace, when `-mask_data` was not given. On an error nothing is
/// undone (NCBI does not either): what was registered and read stays in `db`.
pub(crate) fn process_mask_data(
    db: &mut BuildDatabase<'_>,
    options: &MaskDataOptions<'_>,
    warn: Warn<'_>,
) -> Result<(), MaskDataError> {
    let Some(files) = options.mask_data else {
        return Ok(());
    };
    let file_list = split_commas(truncate_spaces_str(files));
    if file_list.is_empty() {
        return Err(options_error(
            "mask_data option found, but no files were specified.",
        ));
    }

    let mut id_list: Vec<&str> = Vec::new();
    let desc_list: Vec<&str>;

    if let Some(ids) = options.mask_id {
        id_list = split_commas(truncate_spaces_str(ids));
        if file_list.len() != id_list.len() {
            return Err(options_error(
                "the size of mask_id does not match that of mask_data.",
            ));
        }
        // make sure this is not a numeric id
        if id_list.iter().any(|id| is_numeric_int4(id)) {
            return Err(options_error("mask_id can not be numeric."));
        }
    }

    if let Some(descs) = options.mask_desc {
        desc_list = split_commas(truncate_spaces_str(descs));
        if file_list.len() != desc_list.len() {
            return Err(options_error(
                "the size of mask_desc does not match that of mask_data.",
            ));
        }
    } else {
        // description is optional
        desc_list = vec![""; id_list.len()];
    }

    let mut utf8_warned = false;
    for (i, &file) in file_list.iter().enumerate() {
        let path: PathBuf = match options.root {
            Some(root) => root.join(file),
            None => PathBuf::from(file),
        };
        if !std::fs::metadata(&path).is_ok_and(|meta| meta.is_file()) {
            warn(
                format!("Error: [makeblastdb] Ignoring mask file '{file}' as it does not exist.")
                    .as_bytes(),
            );
            continue;
        }

        // The name in the messages is the one given, not the path it was looked up under.
        let mut reader = open_mask_file(&path, Path::new(file))?;

        let mut algo_id: Option<i32> = None;
        loop {
            let next = reader.next_info();
            write_posts(&mut reader, warn, &mut utf8_warned);
            let Some(first_obj) = next? else {
                break;
            };
            let algo_id = match algo_id {
                Some(id) => id,
                None => {
                    let _ = writeln!(db.log(), "Mask file: {file}");
                    let opts = &first_obj.algo_options;
                    let id = if !id_list.is_empty() {
                        db.register_masking_algorithm_id(
                            id_list[i].as_bytes(),
                            desc_list[i].as_bytes(),
                            opts,
                        )?
                    } else {
                        db.register_masking_algorithm(
                            first_obj.algo_program,
                            opts,
                            file.as_bytes(),
                        )?
                    };
                    algo_id = Some(id);
                    id
                }
            };

            let mut masks = first_obj.masks;
            loop {
                // NCBI's list holds the element stated as NULL too (`MaskList::null_element`).
                let listed = !masks.masks.is_empty() || masks.null_element.is_some();
                if db.mask_data_source_mut().is_none() && listed {
                    db.set_mask_data_source(MaskedRangeSet::new());
                }

                for mask in std::mem::take(&mut masks.masks) {
                    let Some(seq_id) = mask.id else {
                        return Err(options_error("Cannot get masked range Seq-id"));
                    };
                    if let Some(set) = db.mask_data_source_mut() {
                        set.insert(algo_id, seq_id, mask.parts)?;
                    }
                }
                // The next mask's `GetId` dereferences an element stated as NULL: refused there.
                if let Some(null) = masks.null_element.take() {
                    return Err(null.refusal().into());
                }

                if !masks.more {
                    break;
                }

                let next = reader.next_list();
                write_posts(&mut reader, warn, &mut utf8_warned);
                masks = next?;
            }
        }
    }
    Ok(())
}

/// Writes on standard error what NCBI's serial reader posted while `reader` read its last object
/// (also when it then failed), in the order posted: before makeblastdb does anything with the
/// object or the error. The warning of a UTF8String tag is `ERR_POST_X_ONCE`, once per run
/// (`utf8_warned`), whichever file it is in.
fn write_posts(reader: &mut MaskFileReader, warn: Warn<'_>, utf8_warned: &mut bool) {
    for post in reader.take_posts() {
        if matches!(post, MaskPost::Utf8Tag(_)) {
            if *utf8_warned {
                continue;
            }
            *utf8_warned = true;
        }
        let line = post.line();
        warn(line.strip_suffix(b"\n").unwrap_or(&line));
    }
}

/// The reader of the mask file at `path`, whose name in messages is `name`. A file that cannot
/// be opened or read is "Unknown encoding" as in NCBI (`MaskFileReader::open_named`).
fn open_mask_file(path: &Path, name: &Path) -> Result<MaskFileReader, MaskDataError> {
    Ok(MaskFileReader::open_named(path, name)?)
}

#[cfg(test)]
#[path = "process_tests.rs"]
mod tests;
