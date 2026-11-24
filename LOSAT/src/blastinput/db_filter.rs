//! The list options of a database search (line DI-8): `-taxids`, `-taxidlist`,
//! `-negative_taxids`, `-negative_taxidlist` and `-no_taxid_expansion` (DI-8a), `-seqidlist`
//! and `-negative_seqidlist` (DI-8b) of BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX. The list
//! is read as NCBI's handler of the database arguments reads it (`s_GetTaxIDList`,
//! `CSeqDBFileGiList`), and the search's database is restricted to it
//! (`seqdb::SeqDb::restrict`).
//!
//! `-gilist` and `-negative_gilist` stay refused; NCBI's exclusions with them are therefore
//! not declared here.

use crate::wording::DATABASE_ERROR;
use std::collections::BTreeSet;
use std::io::Write;

use anyhow::Result;

use crate::blastinput::app::options_error;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::ncbi_str::{string_to_int, ALLOW_LEADING_SPACES, ALLOW_TRAILING_SPACES};
use crate::seqdb::file::DbFile;
use crate::seqdb::lmdb_read::ReadError;
use crate::seqdb::tax_sqlite::{QueryReadFailure, TaxonomySqlite};
use crate::seqdb::{seqidlist, DbList, SearchPath, SeqidList, TaxidList};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2257-2299
// The exclusions of the four lists with one another, with `-subject` and (on the
// programs' `-subject_loc`) with `-subject_loc` are NCBI's argument reader's (`cli::reader`:
// its usage, `Incompatible with argument`, exit status 1).
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2245-2259
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2283-2289
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2302-2318
/// The list options of a search program's command line.
#[derive(Debug, Clone, Default)]
pub struct DbFilterArgs {
    /// `-seqidlist`: a file of Seq-ids (text, or binary from `blastdb_aliastool`).
    pub seqidlist: Option<String>,
    /// `-negative_seqidlist`: a file of Seq-ids to leave out.
    pub negative_seqidlist: Option<String>,
    /// `-taxids`: taxids separated by commas.
    pub taxids: Option<String>,
    /// `-taxidlist`: a file of taxids, one per line.
    pub taxidlist: Option<String>,
    /// `-negative_taxids`: taxids separated by commas.
    pub negative_taxids: Option<String>,
    /// `-negative_taxidlist`: a file of taxids, one per line.
    pub negative_taxidlist: Option<String>,
    /// `-no_taxid_expansion`: the taxids alone, without their descendants.
    pub no_taxid_expansion: bool,
}

impl DbFilterArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, the arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(DbFilterArgs {
            seqidlist: args.given_option("seqidlist")?,
            negative_seqidlist: args.given_option("negative_seqidlist")?,
            taxids: args.given_option("taxids")?,
            taxidlist: args.given_option("taxidlist")?,
            negative_taxids: args.given_option("negative_taxids")?,
            negative_taxidlist: args.given_option("negative_taxidlist")?,
            no_taxid_expansion: args.flag("no_taxid_expansion"),
        })
    }
}

/// The names of the options of `DbFilterArgs`, for the `-subject_loc` exclusions of LOSAT Web's
/// ABI v1, which keeps its own parser (`web_v1_argv`; the other readers are NCBI's).
#[cfg(any(target_arch = "wasm32", test))]
pub const DB_FILTER_OPTIONS: [&str; 7] = [
    "seqidlist",
    "negative_seqidlist",
    "taxids",
    "taxidlist",
    "negative_taxids",
    "negative_taxidlist",
    "no_taxid_expansion",
];

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2426-2431
pub(crate) const NO_TAXONOMY_WARNING: &str =
    "The -taxids command line option requires additional data \
     files. Please see the section 'Taxonomic filtering for BLAST databases' in \
     https://www.ncbi.nlm.nih.gov/books/NBK569839/ for details.\n";

/// NCBI's default name of the taxonomy file (`kDefaultName`, c++/src/objtools/blast/seqdb_reader/tax4blastsqlite.cpp).
pub(crate) const TAXONOMY_FILE: &str = "taxonomy4blast.sqlite3";

impl DbFilterArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2486-2503
    /// Whether the command line gives one of the six lists.
    pub fn given(&self) -> bool {
        self.seqid_list().is_some() || self.list().is_some()
    }

    /// The seqid list option that NCBI reads (before the taxid lists): `(file, is negative)`.
    fn seqid_list(&self) -> Option<(&str, bool)> {
        if let Some(value) = &self.seqidlist {
            Some((value, false))
        } else {
            self.negative_seqidlist
                .as_deref()
                .map(|value| (value, true))
        }
    }

    /// The taxid list option that NCBI reads: `(value, is a file, is negative)`.
    fn list(&self) -> Option<(&str, bool, bool)> {
        if let Some(value) = &self.taxids {
            Some((value, false, false))
        } else if let Some(value) = &self.taxidlist {
            Some((value, true, false))
        } else if let Some(value) = &self.negative_taxids {
            Some((value, false, true))
        } else {
            self.negative_taxidlist
                .as_deref()
                .map(|value| (value, true, true))
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2403-2459
    // Without the taxonomy file NCBI reads every id and keeps none: the set is empty and the
    // search runs on the whole database (DI-8 DECISIONS 2 and 11, the Owner's decision of
    // 2026-10-29). A taxid whose descendants loop back to it (taxid 1, whose parent is 1)
    // never ends NCBI's recursive query; LOSAT visits each taxid once and ends with every
    // taxid below it (DECISIONS 10, `TaxonomySqlite::descendants`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:458-469
    /// The `program-options` parameter of the BLAST archive for the list option of the command
    /// line: NCBI's list object is there for each option, so the parameter is written also when
    /// no id came of it, or the list is a seqid list, without ids. `taxids` are the ids of NCBI's
    /// `set` after the database opened (`SeqDb::list_taxids`), ascending.
    pub fn archive_parameter(
        &self,
        taxids: &[i32],
    ) -> Option<crate::report::archive::options::Param> {
        use crate::report::archive::options::{Param, ParamValue};
        let negative = match (self.seqid_list(), self.list()) {
            (Some((_, negative)), _) | (None, Some((_, _, negative))) => negative,
            (None, None) => return None,
        };
        Some(Param {
            name: if negative {
                "NegativeTaxidList"
            } else {
                "TaxidList"
            },
            value: ParamValue::IntegerList(taxids.to_vec()),
        })
    }

    /// The list of the command line, read (and a taxid list expanded) as NCBI reads it, with
    /// the warnings NCBI writes while reading (a missing taxonomy file, a dropped token of a
    /// seqid list, the taxonomy file's read error) written to `warnings` as the application
    /// `program` (its diagnostic prefix) writes them. `Ok(None)`: no list, or one without an id
    /// (NCBI then searches the whole database).
    pub fn resolve(
        &self,
        settings: &ApplicationSettings,
        program: &str,
        warnings: &mut dyn Write,
    ) -> Result<Option<DbList>> {
        if self.seqid_list().is_none() && self.list().is_none() {
            return Ok(None);
        }
        // The registry's errors are the database's: `open_search_database` reports them.
        let registry = settings.blastdb.clone().unwrap_or(None);
        let path = SearchPath::from_environment(registry.as_deref())
            .map_err(|error| anyhow::anyhow!("{error}"))?;
        if let Some((value, negative)) = self.seqid_list() {
            return self.resolve_seqid_list(value, negative, &path, warnings);
        }
        let Some((value, is_file, negative)) = self.list() else {
            return Ok(None);
        };
        let ids: Vec<Vec<u8>> = if is_file {
            let Some(file) = path.resolve_file(value) else {
                return Err(options_error(&format!("File is not acessible: {value}")));
            };
            let bytes = std::fs::read(&file)
                .map_err(|_| options_error(&format!("File is not acessible: {value}")))?;
            stream_lines(&bytes)
        } else {
            value
                .as_bytes()
                .split(|&byte| byte == b',')
                .filter(|token| !token.is_empty())
                .map(<[u8]>::to_vec)
                .collect()
        };
        let source = if self.no_taxid_expansion {
            None
        } else {
            match path
                .resolve_file(TAXONOMY_FILE)
                .map(|file| DbFile::open(file.as_ref()))
            {
                Some(Ok(file)) => Some(file),
                _ => {
                    warnings.write_all(NO_TAXONOMY_WARNING.as_bytes())?;
                    None
                }
            }
        };
        let taxonomy = match &source {
            Some(file) => match TaxonomySqlite::open(file) {
                Ok(taxonomy) => Some(taxonomy),
                // Every exception of NCBI's constructor, a read that SQLite fails included
                // (`tax_sqlite::QueryReadFailure`; review 1 of line WQA D18b, L3).
                Err(ReadError::Unreadable(_) | ReadError::Format(_)) => {
                    warnings.write_all(NO_TAXONOMY_WARNING.as_bytes())?;
                    None
                }
            },
            None => None,
        };
        let invalid = || options_error("Invalid taxidlist file ");
        let mut tax_ids = BTreeSet::new();
        for id in &ids {
            if crate::ncbi_str::is_blank(id) {
                continue;
            }
            // `TTaxId` is an `int`: `StringToNumeric<TTaxId>` is `NStr::StringToInt`.
            let taxid = string_to_int(id, ALLOW_LEADING_SPACES | ALLOW_TRAILING_SPACES)
                .map_err(|_| invalid())?;
            if self.no_taxid_expansion {
                tax_ids.insert(taxid);
            } else if let Some(taxonomy) = &taxonomy {
                tax_ids.insert(taxid);
                match taxonomy.descendants(taxid) {
                    Ok(descendants) => tax_ids.extend(descendants),
                    // A read that SQLite fails: the statement's destructor posts the error before
                    // the per-id `catch` (`tax_sqlite::QueryReadFailure`; review 1 of line WQA
                    // D18b, L3).
                    Err(ReadError::Unreadable(_)) => {
                        warnings.write_all(&QueryReadFailure::finalize_post(program))?;
                        return Err(invalid());
                    }
                    Err(ReadError::Format(_)) => return Err(invalid()),
                }
            }
        }
        if tax_ids.is_empty() {
            return Ok(None);
        }
        Ok(Some(DbList::Taxids(if negative {
            TaxidList::Negative(tax_ids)
        } else {
            TaxidList::Positive(tax_ids)
        })))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1667-1692
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:254-258
    // The file is looked for on the database search path (`SeqDB_ResolveDbPath`). A file
    // that is not there resolves to the empty name, whose memory map throws a
    // `CFileException` (not a database error): the search program reports it as
    // `Error: <what>` with exit status 255 (`BLAST_UNKNOWN_ERROR`), as is `CMemoryFile`'s exception
    // of a file there that does not open. The reader's own errors (an empty file, a damaged binary
    // list) are `BLAST Database error: ...`, exit status 2.
    /// `CSeqDBFileGiList(fn, eSiList)` of the `-seqidlist` / `-negative_seqidlist` file.
    fn resolve_seqid_list(
        &self,
        value: &str,
        negative: bool,
        path: &SearchPath,
        warnings: &mut dyn Write,
    ) -> Result<Option<DbList>> {
        let bytes = seqidlist::read_resolved(value, path).map_err(|error| {
            match error.ncbi_exception() {
                Some(exception) => anyhow::Error::from(NativeError::exception(exception)),
                None => anyhow::Error::msg(error.to_string()),
            }
        })?;
        let mut cerr = Vec::new();
        let file = seqidlist::read_si_list(&bytes, &mut cerr).map_err(|text| NativeError {
            exit: 2,
            message: format!("{DATABASE_ERROR} {text}\n").into(),
        })?;
        for line in cerr {
            warnings.write_all(line.as_bytes())?;
            warnings.write_all(b"\n")?;
        }
        if file.ids.is_empty() {
            return Ok(None);
        }
        Ok(Some(DbList::Seqids(SeqidList { negative, file })))
    }
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:152-172
// NCBI's reader takes the end-of-line style of the file from its first line (LF, CRLF or a
// lone CR, and both when they mix); LOSAT ends a line at each of them. A carriage return that
// NCBI leaves at an end of a line is white space, which `StringToNumeric` allows around the
// number and `IsBlank` skips, so the taxids are the same. The bytes after the last end are a
// last line when there are any.
/// The lines of a taxid list file.
fn stream_lines(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut at = 0;
    while at < bytes.len() {
        match bytes[at] {
            b'\n' => {
                lines.push(bytes[start..at].to_vec());
                at += 1;
                start = at;
            }
            b'\r' => {
                lines.push(bytes[start..at].to_vec());
                at += if bytes.get(at + 1) == Some(&b'\n') {
                    2
                } else {
                    1
                };
                start = at;
            }
            _ => at += 1,
        }
    }
    if start < bytes.len() {
        lines.push(bytes[start..].to_vec());
    }
    lines
}

#[cfg(test)]
#[path = "db_filter_tests.rs"]
mod tests;
