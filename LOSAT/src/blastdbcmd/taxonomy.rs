//! The taxonomy of `losat blastdbcmd`: the taxids of `-taxids` / `-taxidlist` with their
//! descendants (`x_InitBlastDB_TaxIdList`), whose sequences `retrieve` writes, and
//! `-tax_info`, the taxids of the database (`x_PrintBlastDatabaseTaxInformation`).

use crate::wording::QUERY_OPTIONS_ERROR;
use std::collections::BTreeSet;
use std::io::{BufWriter, Write};

use crate::algorithm::common::database::TaxdbRun;
use crate::blastinput::db_filter::{NO_TAXONOMY_WARNING, TAXONOMY_FILE};
use crate::ncbi_app::args::{arg_message, Args};
use crate::ncbi_app::files::{
    flush_before_diagnostic, input_lines, read_input_file, Cin, LineEnds, Ostream,
};
use crate::ncbi_str::{string_to_int, ALLOW_LEADING_SPACES, ALLOW_TRAILING_SPACES};
use crate::seqdb::file::DbFile;
use crate::seqdb::lmdb_read::ReadError;
use crate::seqdb::tax_sqlite::{QueryReadFailure, TaxonomySqlite};
use crate::seqdb::taxdb::{TaxLookup, TaxNames};
use crate::seqdb::{SearchPath, SeqDb, SeqDbError};

use super::{database_error, refusal, Diag};

/// What ends the reading of the taxids before the database opens: the exit status, the
/// diagnostics already posted.
pub(super) struct Ended(pub i32);

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:319-357
// `NStr::Split` keeps the empty pieces of `1,,2` and `1,` and gives none for the empty
// string (probes J101, J102, S36 of line BDC). The lines of the file are those of
// `NcbiGetlineEOL` (`args::input_lines`), as `-entry_batch` reads them: on Linux a line
// ends at `\n` alone, and a `\r` stays and reads as white space around the number (J109).
// Unlike a search program's list (`db_filter`), the id is kept without the taxonomy file
// (T02), and a line of blanks is read as a number. An id that is not a number escapes
// NCBI's `Run` as a `CStringException` with the path of NCBI's build in its text (exit
// status 255: J101, J102, J107, J108, J131, S35), which LOSAT refuses (line BDC DECISIONS
// 8). `-taxids 1` never ends in NCBI (the parent of taxid 1 is 1, probe J120); LOSAT's walk
// visits each taxid once and ends with every taxid below it (`TaxonomySqlite::descendants`,
// DI-8 DECISIONS 10).
/// The taxids of `-taxids` / `-taxidlist` (`-` reads `stdin`), with their descendants in the
/// taxonomy file of the search path unless `-no_taxid_expansion`; NCBI's warning when that file
/// is missing or not one NCBI reads goes to `diag`.
pub(super) fn tax_id_list(
    args: &Args,
    path: &SearchPath,
    diag: &mut Diag,
    stdin: &mut Cin<'_>,
) -> Result<BTreeSet<i32>, Ended> {
    let ids: Vec<Vec<u8>> = match (args.value("taxids"), args.value("taxidlist")) {
        (Some(""), _) => Vec::new(),
        (Some(input), _) => input.split(',').map(|id| id.as_bytes().to_vec()).collect(),
        (None, Some(file)) => {
            let Ok(bytes) = read_input_file("taxidlist", file, stdin) else {
                // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
                diag.post(&format!(
                    "Command line argument error: {}",
                    arg_message("taxidlist", "File is not accessible", file)
                ));
                return Err(Ended(1));
            };
            input_lines(&bytes, LineEnds::HOST)
                .into_iter()
                .map(<[u8]>::to_vec)
                .collect()
        }
        (None, None) => unreachable!("the caller reads a taxid list"),
    };
    let source = if args.given("no_taxid_expansion") {
        None
    } else {
        match path
            .resolve_file(TAXONOMY_FILE)
            .map(|file| DbFile::open(file.as_ref()))
        {
            Some(Ok(file)) => Some(file),
            _ => {
                diag.post(NO_TAXONOMY_WARNING.trim_end_matches('\n'));
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
                diag.post(NO_TAXONOMY_WARNING.trim_end_matches('\n'));
                None
            }
        },
        None => None,
    };
    let mut tax_ids = BTreeSet::new();
    for id in &ids {
        // `TTaxId` is an `int`: `StringToNumeric<TTaxId>` is `NStr::StringToInt`.
        let Ok(taxid) = string_to_int(id, ALLOW_LEADING_SPACES | ALLOW_TRAILING_SPACES) else {
            return Err(Ended(diag.refuse(&refusal(&format!(
                "the taxid '{}' is not a number: NCBI ends with an exception that names its \
                 source files (exit status 255)",
                String::from_utf8_lossy(id)
            )))));
        };
        tax_ids.insert(taxid);
        if let Some(taxonomy) = &taxonomy {
            let descendants = taxonomy.descendants(taxid).map_err(|error| match error {
                // A read that SQLite fails: the statement's destructor posts the error, then the
                // step's exception ends the run (`tax_sqlite::QueryReadFailure`; review 1 of line
                // WQA D18b, L3).
                ReadError::Unreadable(_) => {
                    diag.write_text(&QueryReadFailure::finalize_post("blastdbcmd"));
                    Ended(database_error(
                        diag,
                        SeqDbError::Ncbi(QueryReadFailure::step()),
                    ))
                }
                ReadError::Format(reason) => Ended(diag.refuse(&refusal(&format!(
                    "{TAXONOMY_FILE}: the descendants of taxid {taxid} are not readable \
                     ({reason})"
                )))),
            })?;
            tax_ids.extend(descendants);
        }
    }
    Ok(tax_ids)
}

/// A field of `-tax_info` (`CPrintTaxFields`' enum).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Field {
    TaxId,
    SciName,
    CommonName,
    SuperKingdom,
    BlastName,
    NumSeqs,
}

/// `CPrintTaxFields`: the fields of the `-outfmt` of `-tax_info` and the text before each,
/// with the text after the last.
#[derive(Debug, PartialEq, Eq)]
struct TaxFields {
    fields: Vec<Field>,
    /// `m_Seperators`: one more than `fields` for a format without `%f`; `%f` pushes its
    /// own, so `separators[i]` precedes field `i` and the last one ends the line.
    separators: Vec<Vec<u8>>,
    need_names: bool,
    need_num_seqs: bool,
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:592-653
// Ported as it reads, the quirks with it (probes T10-T12, J140-J159 of line BDC): `%f`
// pushes the text before it and a tab per field, then the text after the `switch` again,
// so after `%f` the fields shift onto the tabs (`%f %n` is `T\tS\tL\tK\tB\tN`); `%%` adds one
// `%` and the second `%` starts a specifier (`%%T` is `%` and the taxid, `%%x` the text
// `%%x`); a `%` that ends the format reads the string's NUL and adds `%` and a NUL byte;
// any other letter adds itself and its `%` as text.
impl TaxFields {
    /// `CPrintTaxFields` with `out` and `fmt`; `Err` for NCBI's
    /// `Invalid format options for tax_info.`
    fn parse(fmt: &[u8]) -> Result<Self, ()> {
        let at = |i: usize| fmt.get(i).copied().unwrap_or(0);
        let mut fields = Vec::new();
        let mut separators = Vec::new();
        let (mut need_names, mut need_num_seqs) = (false, false);
        let mut sp = Vec::new();
        let mut i = 0;
        while i < fmt.len() {
            if fmt[i] != b'%' {
                sp.push(fmt[i]);
                i += 1;
                continue;
            }
            if at(i + 1) == b'%' {
                sp.push(b'%');
                i += 1;
                continue;
            }
            i += 1;
            match at(i) {
                b'f' => {
                    need_names = true;
                    separators.push(sp.clone());
                    for field in [
                        Field::TaxId,
                        Field::SciName,
                        Field::CommonName,
                        Field::SuperKingdom,
                        Field::BlastName,
                    ] {
                        fields.push(field);
                        separators.push(b"\t".to_vec());
                    }
                }
                b'T' => fields.push(Field::TaxId),
                b'S' => {
                    fields.push(Field::SciName);
                    need_names = true;
                }
                b'L' => {
                    fields.push(Field::CommonName);
                    need_names = true;
                }
                b'K' => {
                    fields.push(Field::SuperKingdom);
                    need_names = true;
                }
                b'B' => {
                    fields.push(Field::BlastName);
                    need_names = true;
                }
                b'n' => {
                    fields.push(Field::NumSeqs);
                    need_num_seqs = true;
                }
                other => {
                    sp.push(b'%');
                    sp.push(other);
                    i += 1;
                    continue;
                }
            }
            separators.push(std::mem::take(&mut sp));
            i += 1;
        }
        separators.push(sp);
        if fields.is_empty() {
            return Err(());
        }
        Ok(TaxFields {
            fields,
            separators,
            need_names,
            need_num_seqs,
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:657-697
    /// `PrintEntry`: the line of `taxid` with its `names` and its number of sequences.
    fn line(&self, taxid: i32, names: &LineNames, num_seqs: usize, out: &mut Vec<u8>) {
        let name = |pick: fn(&TaxNames) -> &[u8]| -> &[u8] {
            match names {
                LineNames::Named(names) => pick(names),
                LineNames::Empty => b"",
                LineNames::NotAvailable => b"N/A",
            }
        };
        for (i, field) in self.fields.iter().enumerate() {
            out.extend_from_slice(&self.separators[i]);
            match field {
                Field::TaxId => out.extend_from_slice(taxid.to_string().as_bytes()),
                Field::SciName => out.extend_from_slice(name(|n| &n.scientific_name)),
                Field::CommonName => out.extend_from_slice(name(|n| &n.common_name)),
                // `tax_info` refuses a kingdom past NCBI's mapping of taxdb.btd.
                Field::SuperKingdom => {
                    out.extend_from_slice(name(|n| n.super_kingdom.as_deref().unwrap_or_default()))
                }
                Field::BlastName => out.extend_from_slice(name(|n| &n.blast_name)),
                Field::NumSeqs => out.extend_from_slice(num_seqs.to_string().as_bytes()),
            }
        }
        out.extend_from_slice(self.separators.last().map_or(&[][..], Vec::as_slice));
        out.push(b'\n');
    }
}

/// The names of a `-tax_info` line (`SSeqDBTaxInfo` after `GetTaxNames`).
enum LineNames {
    /// The lookup found them.
    Named(TaxNames),
    /// The lookup found the taxid and no names (`info.taxid` set, the names empty).
    Empty,
    /// The lookup did not find the taxid, or the taxid is 0: NCBI's `N/A`.
    NotAvailable,
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:699-733
// NCBI's exceptions escape to `Run`: the format's (`CInputException`, exit status 1), the
// reader's (`CSeqDBException`, exit status 2) after the lines already written. Taxid 0 (the
// sequences without one) is a taxid like any other; taxdb has no names for it (T12). The
// `info.taxid` that `GetTaxNames` sets when its search finds the taxid (seqdbtax.cpp:258-259)
// stays set where it then returns false (the offset error, a record not of four fields): such a
// line has the taxid and empty names, `N/A` only where the search did not find it, and always
// for taxid 0 (probe of line RG T5b, `T/t5b/out/extra`). The lookups' texts (`ERR_POST`: of
// opening the database, the offset error) come before the line, after those written (the
// `endl` of the count, `cerr`'s tie).
/// `-tax_info` of `db`: the number of its taxids, then a line for each, ascending, to `out`;
/// the exit status.
pub(super) fn tax_info(db: &SeqDb, args: &Args, diag: &mut Diag, out: &mut Ostream<'_>) -> i32 {
    let fmt = args.value("outfmt").unwrap_or("%f");
    let Ok(fields) = TaxFields::parse(fmt.as_bytes()) else {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-176
        diag.post(&format!(
            "{QUERY_OPTIONS_ERROR} Invalid format options for tax_info."
        ));
        // Exception (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, exception 12).
        diag.post(crate::wording::USER_MANUAL_FOOTER);
        return 1;
    };
    let lookup = match db.tax_id_lookup() {
        Ok(lookup) => lookup,
        Err(error) => return database_error(diag, error),
    };
    let tax_ids = match lookup.db_tax_ids() {
        Ok(tax_ids) => tax_ids,
        Err(error) => return database_error(diag, error),
    };
    let taxdb = TaxdbRun::sequential();
    let need_kingdom = fields.fields.contains(&Field::SuperKingdom);
    let mut sink = BufWriter::new(out);
    let _ = writeln!(sink, "# of Tax IDs in Database: {}", tax_ids.len());
    let mut line = Vec::new();
    for &taxid in &tax_ids {
        let names = if fields.need_names {
            let mut messages = Vec::new();
            let lookup = taxdb.lookup(db, taxid, &mut messages);
            for message in &messages {
                if flush_before_diagnostic(&mut sink) {
                    return 0;
                }
                diag.err_post(message);
            }
            // A lookup whose result NCBI reads past its mapping of taxdb.btd is refused where
            // the line writes it: names of no defined result, a kingdom past the mapping. Taxid
            // 0's line writes `N/A` whatever its lookup gives (NCBI's `ZERO_TAX_ID` test above).
            let lookup = lookup.and_then(|lookup| match lookup {
                lookup if taxid == 0 => Ok(lookup),
                TaxLookup::PastMapping(past) => Err(past.error("names")),
                TaxLookup::Named(names) if need_kingdom => match &names.super_kingdom {
                    Err(past) => Err(past.error("super kingdom")),
                    Ok(_) => Ok(TaxLookup::Named(names)),
                },
                lookup => Ok(lookup),
            });
            match lookup {
                Ok(_) if taxid == 0 => LineNames::NotAvailable,
                Ok(TaxLookup::Named(names)) => LineNames::Named(names),
                Ok(TaxLookup::Unnamed) => LineNames::Empty,
                Ok(TaxLookup::NotInIndex | TaxLookup::PastMapping(_)) => LineNames::NotAvailable,
                Err(error) => {
                    if flush_before_diagnostic(&mut sink) {
                        return 0;
                    }
                    return database_error(diag, error);
                }
            }
        } else {
            LineNames::NotAvailable
        };
        let num_seqs = if fields.need_num_seqs {
            match lookup.tax_ids_to_oids(&BTreeSet::from([taxid])) {
                Ok(oids) => oids.len(),
                Err(error) => {
                    if flush_before_diagnostic(&mut sink) {
                        return 0;
                    }
                    return database_error(diag, error);
                }
            }
        } else {
            0
        };
        line.clear();
        fields.line(taxid, &names, num_seqs, &mut line);
        // `out` is the run's `ncbi_app::files::Ostream`: a failed write does nothing, as in
        // NCBI's stream (no exception is enabled on it); a closed pipe ends the program, where
        // SIGPIPE ends NCBI's.
        if sink.write_all(&line).is_err() {
            return 0;
        }
    }
    let _ = sink.flush();
    0
}

#[cfg(test)]
#[path = "taxonomy_tests.rs"]
mod tests;
