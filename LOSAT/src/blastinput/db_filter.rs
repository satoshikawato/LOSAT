//! The list options of a database search: `-taxids`, `-taxidlist`, `-negative_taxids`,
//! `-negative_taxidlist` and `-no_taxid_expansion` (DI-8a), `-seqidlist` and `-negative_seqidlist`
//! (DI-8b), `-gilist` and `-negative_gilist` across all six database search programs, and
//! `-ipglist` and `-negative_ipglist` on protein-database search programs (`blastp`, `blastx`,
//! `psiblast`). The list is read as NCBI's handler of the database arguments reads it
//! (`s_GetTaxIDList`, `CSeqDBFileGiList`), and the search's database is restricted to it
//! (`seqdb::SeqDb::restrict`).

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
use crate::seqdb::{id_list, seqidlist, DbList, GiList, PigList, SearchPath, SeqidList, TaxidList};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2257-2326
// The exclusions among the list options, with `-subject` and (on the programs' `-subject_loc`)
// with `-subject_loc` are NCBI's argument reader's (`cli::reader`: its usage,
// `Incompatible with argument`, exit status 1). Note that `-ipglist` and `-negative_ipglist` are
// not in NCBI's `kBlastDBFilteringOptions`, so they can be passed alongside another list option
// and `ExtractAlgorithmOptions`'s branch order decides which one takes effect.
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2211-2259
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2282-2289
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2302-2326
/// The list options of a search program's command line.
#[derive(Debug, Clone, Default)]
pub struct DbFilterArgs {
    /// `-gilist`: a file of GIs (text or binary).
    pub gilist: Option<String>,
    /// `-negative_gilist`: a file of GIs to leave out.
    pub negative_gilist: Option<String>,
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
    /// `-ipglist`: a file of IPGs/PIGs (text or binary; protein databases only).
    pub ipglist: Option<String>,
    /// `-negative_ipglist`: a file of IPGs/PIGs to leave out (protein databases only).
    pub negative_ipglist: Option<String>,
    /// `-no_taxid_expansion`: the taxids alone, without their descendants.
    pub no_taxid_expansion: bool,
}

impl DbFilterArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args` for a nucleotide-database search program (`blastn`, `tblastn`,
    /// `tblastx`), which has no `-ipglist` or `-negative_ipglist`.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Self::from_args_inner(args, false)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:63
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:64
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:59
    /// The options of `args` for a protein-database search program (`blastp`, `blastx`,
    /// `psiblast`), which enables `-ipglist` and `-negative_ipglist` via
    /// `SetIPGFilteringSupport`.
    pub(crate) fn from_args_with_ipg(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Self::from_args_inner(args, true)
    }

    fn from_args_inner(
        args: &crate::ncbi_app::args::Args,
        with_ipg: bool,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(DbFilterArgs {
            gilist: args.given_option("gilist")?,
            negative_gilist: args.given_option("negative_gilist")?,
            seqidlist: args.given_option("seqidlist")?,
            negative_seqidlist: args.given_option("negative_seqidlist")?,
            taxids: args.given_option("taxids")?,
            taxidlist: args.given_option("taxidlist")?,
            negative_taxids: args.given_option("negative_taxids")?,
            negative_taxidlist: args.given_option("negative_taxidlist")?,
            ipglist: if with_ipg {
                args.given_option("ipglist")?
            } else {
                None
            },
            negative_ipglist: if with_ipg {
                args.given_option("negative_ipglist")?
            } else {
                None
            },
            no_taxid_expansion: args.flag("no_taxid_expansion"),
        })
    }
}

/// The frozen set of `DbFilterArgs` options from LOSAT Web's ABI v1 (which keeps its own parser
/// `web_v1_argv` and does not accept `-gilist`, `-negative_gilist`, `-ipglist`, or
/// `-negative_ipglist`), used for v1 `-subject_loc` exclusions.
#[cfg(any(losat_web_abi_v1, test))]
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

/// Which list option of `DbFilterArgs` wins NCBI's `ExtractAlgorithmOptions` branch chain.
#[derive(Clone, Copy)]
enum SelectedList<'a> {
    Gi {
        value: &'a str,
        negative: bool,
    },
    Seqid {
        value: &'a str,
        negative: bool,
    },
    Taxid {
        value: &'a str,
        is_file: bool,
        negative: bool,
    },
    Pig {
        value: &'a str,
        negative: bool,
    },
}

impl DbFilterArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2478-2511
    /// The first given list option in NCBI's `ExtractAlgorithmOptions` order.
    fn selected(&self) -> Option<SelectedList<'_>> {
        if let Some(value) = &self.gilist {
            Some(SelectedList::Gi {
                value,
                negative: false,
            })
        } else if let Some(value) = &self.negative_gilist {
            Some(SelectedList::Gi {
                value,
                negative: true,
            })
        } else if let Some(value) = &self.seqidlist {
            Some(SelectedList::Seqid {
                value,
                negative: false,
            })
        } else if let Some(value) = &self.negative_seqidlist {
            Some(SelectedList::Seqid {
                value,
                negative: true,
            })
        } else if let Some(value) = &self.taxids {
            Some(SelectedList::Taxid {
                value,
                is_file: false,
                negative: false,
            })
        } else if let Some(value) = &self.taxidlist {
            Some(SelectedList::Taxid {
                value,
                is_file: true,
                negative: false,
            })
        } else if let Some(value) = &self.negative_taxids {
            Some(SelectedList::Taxid {
                value,
                is_file: false,
                negative: true,
            })
        } else if let Some(value) = &self.negative_taxidlist {
            Some(SelectedList::Taxid {
                value,
                is_file: true,
                negative: true,
            })
        } else if let Some(value) = &self.ipglist {
            Some(SelectedList::Pig {
                value,
                negative: false,
            })
        } else {
            self.negative_ipglist
                .as_deref()
                .map(|value| SelectedList::Pig {
                    value,
                    negative: true,
                })
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2478-2511
    /// Whether the command line gives one of the database list options.
    pub fn given(&self) -> bool {
        self.selected().is_some()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2403-2459
    // Without the taxonomy file NCBI reads every id and keeps none: the set is empty and the
    // search runs on the whole database (DI-8 DECISIONS 2 and 11, the Owner's decision of
    // 2026-10-29). A taxid whose descendants loop back to it (taxid 1, whose parent is 1)
    // never ends NCBI's recursive query; LOSAT visits each taxid once and ends with every
    // taxid below it (DECISIONS 10, `TaxonomySqlite::descendants`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:444-498
    /// The `program-options` parameters of the BLAST archive and search strategy for the list
    /// option of the command line: for a GI list with at least one GI, `GiList` or
    /// `NegativeGiList` (`BigIntegerList`) followed by an empty `TaxidList` or
    /// `NegativeTaxidList`; for every other list (including a 0-GI file, seqid, taxid, or IPG
    /// list), `TaxidList` or `NegativeTaxidList` alone. `gis` are the GIs of NCBI's list object
    /// after the database opened (`SeqDb::list_gis`), `taxids` are its taxids
    /// (`SeqDb::list_taxids`).
    pub fn archive_parameters(
        &self,
        gis: &[i64],
        taxids: &[i32],
    ) -> Vec<crate::report::archive::options::Param> {
        use crate::report::archive::options::{Param, ParamValue};
        let Some(selected) = self.selected() else {
            return Vec::new();
        };
        let (is_gi, negative) = match selected {
            SelectedList::Gi { negative, .. } => (true, negative),
            SelectedList::Seqid { negative, .. }
            | SelectedList::Taxid { negative, .. }
            | SelectedList::Pig { negative, .. } => (false, negative),
        };
        let mut params = Vec::with_capacity(2);
        if is_gi && !gis.is_empty() {
            params.push(Param {
                name: if negative { "NegativeGiList" } else { "GiList" },
                value: ParamValue::BigIntegerList(gis.to_vec()),
            });
        }
        params.push(Param {
            name: if negative {
                "NegativeTaxidList"
            } else {
                "TaxidList"
            },
            value: ParamValue::IntegerList(taxids.to_vec()),
        });
        params
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
        let Some(selected) = self.selected() else {
            return Ok(None);
        };
        let path = SearchPath::from_environment(settings.blastdb.as_deref());
        let (value, is_file, negative) = match selected {
            SelectedList::Gi { value, negative } => {
                return self.resolve_gi_list(value, negative, &path);
            }
            SelectedList::Seqid { value, negative } => {
                return self.resolve_seqid_list(value, negative, &path, warnings);
            }
            SelectedList::Pig { value, negative } => {
                return self.resolve_pig_list(value, negative, &path);
            }
            SelectedList::Taxid {
                value,
                is_file,
                negative,
            } => (value, is_file, negative),
        };
        let ids: Vec<Vec<u8>> = if is_file {
            let Some(file) = path.resolve_file(value.as_bytes()) else {
                return Err(options_error(&format!("File is not acessible: {value}")));
            };
            let bytes = std::fs::read(crate::utils::os_bytes::path_of_bytes(&file))
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
                .resolve_file(TAXONOMY_FILE.as_bytes())
                .map(|file| DbFile::open(&crate::utils::os_bytes::path_of_bytes(&file)))
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:254-257,265-276,290-291
    // The file is looked for on the database search path (`SeqDB_ResolveDbPath`). A file
    // that is not there resolves to the empty name, whose memory map throws a
    // `CFileException` (not a database error): the search program reports it as
    // `Error: <what>` with exit status 255 (`BLAST_UNKNOWN_ERROR`), as is `CMemoryFile`'s exception
    // of a file there that does not open. The reader's own errors (an empty file, a damaged binary
    // list) are `BLAST Database error: ...`, exit status 2.
    /// Reads a resolved list file (`SeqDB_ResolveDbPath`) with NCBI's `CFileException` mapping.
    fn read_list_file(value: &str, path: &SearchPath) -> Result<Vec<u8>> {
        seqidlist::read_resolved(value.as_bytes(), path).map_err(|error| {
            match error.ncbi_exception() {
                Some(exception) => anyhow::Error::from(NativeError::exception(exception)),
                None => anyhow::Error::msg(error.to_string()),
            }
        })
    }

    /// `CSeqDBFileGiList(fn, eGiList)` of the `-gilist` / `-negative_gilist` file.
    fn resolve_gi_list(
        &self,
        value: &str,
        negative: bool,
        path: &SearchPath,
    ) -> Result<Option<DbList>> {
        let bytes = Self::read_list_file(value, path)?;
        let gis = id_list::read_gi_list(&bytes).map_err(|text| NativeError {
            exit: 2,
            message: format!("{DATABASE_ERROR} {text}\n").into(),
        })?;
        if gis.is_empty() {
            return Ok(None);
        }
        Ok(Some(DbList::Gis(GiList { negative, gis })))
    }

    /// `CSeqDBFileGiList(fn, ePigList)` of the `-ipglist` / `-negative_ipglist` file.
    fn resolve_pig_list(
        &self,
        value: &str,
        negative: bool,
        path: &SearchPath,
    ) -> Result<Option<DbList>> {
        let bytes = Self::read_list_file(value, path)?;
        let pigs = id_list::read_pig_list(&bytes).map_err(|text| NativeError {
            exit: 2,
            message: format!("{DATABASE_ERROR} {text}\n").into(),
        })?;
        if pigs.is_empty() {
            return Ok(None);
        }
        Ok(Some(DbList::Pigs(PigList { negative, pigs })))
    }

    /// `CSeqDBFileGiList(fn, eSiList)` of the `-seqidlist` / `-negative_seqidlist` file.
    fn resolve_seqid_list(
        &self,
        value: &str,
        negative: bool,
        path: &SearchPath,
        warnings: &mut dyn Write,
    ) -> Result<Option<DbList>> {
        let bytes = Self::read_list_file(value, path)?;
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
