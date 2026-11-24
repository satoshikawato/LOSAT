//! blastdbcheck's test of an alias file (`CAliasTest`): its lines, the files that its lists name,
//! the pairing of its OID range and lists, and the database that it makes.

use crate::ncbi_str::{is_blank, string_to_int, string_to_int8};
use crate::seqdb::{SearchPath, SeqDb, SeqType};
use crate::utils::os_bytes::path_from_bytes;

use super::log::{Log, Stop, BRIEF, DETAILS};
use super::path::{make_path, FileParts};
use super::work::{caught, unknown_text};
use super::PROGRAM;

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-182
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-298
// `CStreamLineReader` of the default end-of-line style (`eEOL_unknown`) learns the style from the
// first line and switches as it meets the other byte: whichever of CR and LF the style reads up
// to, a line that holds the other is cut there and its rest read again. So each CR and each LF
// ends a line; a CR LF ends one, or ends one and an empty line, and the test skips empty lines
// (as it skips the empty piece after a last line end here). A file whose bytes cannot be read
// gives no line (the stream's state ends the loop).
/// The lines of an alias file as blastdbcheck reads them: its bytes cut at each CR and LF.
fn lines(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    bytes.split(|&byte| byte == b'\n' || byte == b'\r')
}

/// The values of an alias file's keywords that the test checks after its lines; -1 for one that
/// is absent or not a number, 0 for a list that is present.
struct Values {
    length: i64,
    oidlist: i32,
    gilist: i32,
    nseq: i32,
    first_oid: i32,
    last_oid: i32,
    maxoid: i32,
}

/// The value of an `int` expression of an alias file's values, exact and as the machine's
/// arithmetic wraps it.
#[derive(Clone, Copy)]
struct Both {
    exact: i128,
    wrapped: i64,
}

impl Both {
    fn int(value: i32) -> Self {
        Both {
            exact: value.into(),
            wrapped: value.into(),
        }
    }

    fn int8(value: i64) -> Self {
        Both {
            exact: value.into(),
            wrapped: value,
        }
    }

    /// `self + other` in `int` (`int8`: in `Int8`).
    fn add(self, other: Self, int8: bool) -> Self {
        Both {
            exact: self.exact + other.exact,
            wrapped: wrap(self.wrapped.wrapping_add(other.wrapped), int8),
        }
    }

    /// `self - other` in `int`.
    fn sub(self, other: Self) -> Self {
        Both {
            exact: self.exact - other.exact,
            wrapped: wrap(self.wrapped.wrapping_sub(other.wrapped), false),
        }
    }

    /// `self * other` in `int` (`int8`: in `Int8`).
    fn mul(self, other: Self, int8: bool) -> Self {
        Both {
            exact: self.exact * other.exact,
            wrapped: wrap(self.wrapped.wrapping_mul(other.wrapped), int8),
        }
    }
}

/// `value` cut to `int` (unless `int8`).
fn wrap(value: i64, int8: bool) -> i64 {
    if int8 {
        value
    } else {
        i64::from(value as i32)
    }
}

// The comparisons NCBI makes in `int` (and in `Int8` for OIDLIST) overflow for large values,
// whose result C leaves undefined: NCBI's build wraps them (row RG-D2b a_ov1: FIRST_OID 65535
// and LAST_OID 65535 multiply to 0 in `int`). Whether LOSAT gives the exact result there or the
// wrapped one is the Owner's question Q-D2b-OV; until then LOSAT refuses an alias file whose
// result differs between the two, and gives the result where it does not (a_ov3).
/// The truth of `test` on an expression's value, where the exact and the wrapped values agree;
/// LOSAT's refusal where they do not (`what` names the check).
fn decide(value: Both, test: impl Fn(i128) -> bool, name: &[u8], what: &str) -> Result<bool, Stop> {
    let exact = test(value.exact);
    if exact == test(value.wrapped.into()) {
        return Ok(exact);
    }
    Err(Stop::Refused(format!(
        "the alias file {} whose values overflow NCBI's arithmetic in its check of {what} (the \
         result of the overflow, which C leaves undefined, decides the check), which is not \
         supported by LOSAT's {PROGRAM}",
        String::from_utf8_lossy(name)
    )))
}

/// The test of the alias file `name`: its faults.
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:698-807
pub(super) fn alias_file_test(
    log: &mut Log<'_>,
    name: &str,
    path: &SearchPath,
) -> Result<u32, Stop> {
    let name = name.as_bytes();
    let file = FileParts::of(name);
    let prefix: Vec<u8> = [b"  ", name, b" / AliasFileTest: "].concat();
    let mut faults = 0u32;
    match file.file_length() {
        None => {
            log.line(BRIEF, &[&prefix, b"  [ERROR] alias file does not exist"])?;
            return Ok(1);
        }
        Some(0) => {
            log.line(BRIEF, &[&prefix, b"  [ERROR] alias file has zero length"])?;
            return Ok(1);
        }
        Some(_) => {}
    }
    let mut values = Values {
        length: -1,
        oidlist: -1,
        gilist: -1,
        nseq: -1,
        first_oid: -1,
        last_oid: -1,
        maxoid: -1,
    };
    let bytes = std::fs::read(path_from_bytes(file.path)).unwrap_or_default();
    for line in lines(&bytes) {
        if line.first() == Some(&b'#') || is_blank(line) {
            continue;
        }
        let tokens: Vec<&[u8]> = line
            .split(|&byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r'))
            .filter(|token| !token.is_empty())
            .collect();
        // A line that is not blank holds a byte that is not a separator.
        let Some(&keyword) = tokens.first() else {
            continue;
        };
        if tokens.len() == 1 {
            log.line(
                BRIEF,
                &[
                    &prefix,
                    b"  [ERROR] no value(s) found for keyword ",
                    keyword,
                ],
            )?;
            faults += 1;
            continue;
        }
        match keyword {
            b"DBLIST" | b"TITLE" | b"MASKLIST" | b"STATS_TOTLEN" | b"STATS_NSEQ" => {}
            b"TAXIDLIST" => faults += check_files(log, &prefix, &file.dir, &tokens)?,
            b"OIDLIST" => {
                values.oidlist = 0;
                faults += check_files(log, &prefix, &file.dir, &tokens)?;
            }
            b"GILIST" => {
                values.gilist = 0;
                faults += check_files(log, &prefix, &file.dir, &tokens)?;
            }
            b"NSEQ" | b"MEMB_BIT" | b"FIRST_OID" | b"LAST_OID" | b"MAXOID" | b"MAXLEN" => {
                let value = match string_to_int(tokens[1], 0) {
                    Ok(value) => value,
                    Err(_) => {
                        not_a_number(log, &prefix, keyword)?;
                        faults += 1;
                        -1
                    }
                };
                match keyword {
                    b"NSEQ" => values.nseq = value,
                    b"FIRST_OID" => values.first_oid = value,
                    b"LAST_OID" => values.last_oid = value,
                    b"MAXOID" => values.maxoid = value,
                    _ => {}
                }
            }
            b"LENGTH" => {
                values.length = match string_to_int8(tokens[1], 0) {
                    Ok(value) => value,
                    Err(_) => {
                        not_a_number(log, &prefix, keyword)?;
                        faults += 1;
                        -1
                    }
                };
            }
            _ => {
                log.line(
                    BRIEF,
                    &[&prefix, b"  [ERROR] unknown keyword encountered: ", keyword],
                )?;
                faults += 1;
            }
        }
    }
    faults += check_values(log, &prefix, name, &values)?;
    faults += open(log, &prefix, &file, path)?;
    Ok(faults)
}

/// The line of a value that is not a number for `keyword`.
fn not_a_number(log: &mut Log<'_>, prefix: &[u8], keyword: &[u8]) -> Result<(), Stop> {
    log.line(
        BRIEF,
        &[
            prefix,
            b"  [ERROR] could not convert value to number for ",
            keyword,
        ],
    )
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:832-871
/// The check of the files that a list keyword's values name in the directory `dir` of the alias
/// file (`x_CheckFile`): its faults.
fn check_files(
    log: &mut Log<'_>,
    prefix: &[u8],
    dir: &[u8],
    tokens: &[&[u8]],
) -> Result<u32, Stop> {
    let keyword = tokens[0];
    let mut faults = 0;
    for &token in &tokens[1..] {
        let joined = make_path(dir, token);
        let file = FileParts::of(&joined);
        let warning: &[u8] = match keyword {
            b"OIDLIST" if file.ext != b".msk" => b"  [WARNING] oidlist file ",
            b"GILIST" if file.ext != b".gil" => b"  [WARNING] gilist file ",
            _ => b"",
        };
        if !warning.is_empty() {
            let extension: &[u8] = if keyword == b"OIDLIST" {
                b".msk"
            } else {
                b".gil"
            };
            log.line(
                DETAILS,
                &[
                    prefix,
                    warning,
                    token,
                    b" does not have ",
                    extension,
                    b" extension",
                ],
            )?;
        }
        let problem: &[u8] = match file.file_length() {
            None => b" does not exist",
            Some(0) => b" has zero length",
            Some(_) => continue,
        };
        log.line(
            BRIEF,
            &[
                prefix,
                b"  [ERROR] file ",
                token,
                b" referenced in keyword ",
                keyword,
                problem,
            ],
        )?;
        faults += 1;
    }
    Ok(faults)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:771-792
/// The checks of the OID range, the lists and their counts after the lines: its faults.
fn check_values(
    log: &mut Log<'_>,
    prefix: &[u8],
    name: &[u8],
    values: &Values,
) -> Result<u32, Stop> {
    let one = Both::int(1);
    let first = Both::int(values.first_oid);
    let last = Both::int(values.last_oid);
    let nseq = Both::int(values.nseq);
    let mut faults = 0;
    let unpaired = first.add(one, false).mul(last.add(one, false), false);
    if decide(unpaired, |v| v == 0, name, "FIRST_OID and LAST_OID")? {
        let sum = first.add(last, false);
        if decide(sum, |v| v != -2, name, "FIRST_OID and LAST_OID")? {
            log.line(
                BRIEF,
                &[prefix, b"  [ERROR] FIRST_OID not paired with LAST_OID"],
            )?;
            faults += 1;
        }
    } else {
        let span = nseq.add(first, false).sub(last).sub(one);
        let inconsistent = nseq.add(one, false).mul(span, false);
        if decide(
            inconsistent,
            |v| v != 0,
            name,
            "NSEQ with FIRST_OID and LAST_OID",
        )? {
            log.line(
                BRIEF,
                &[
                    prefix,
                    b"  [ERROR] (FIRST_OID, LAST_OID) is not consistent with NSEQ",
                ],
            )?;
            faults += 1;
        }
    }
    if values.oidlist + values.gilist >= 0 {
        log.line(
            BRIEF,
            &[
                prefix,
                b"  [ERROR] OIDLIST and GILIST cannot be present in the same alias file",
            ],
        )?;
        faults += 1;
    }
    if values.oidlist != -1 {
        let counts = Both::int(values.maxoid)
            .add(one, false)
            .mul(Both::int8(values.length).add(Both::int8(1), true), true)
            .mul(nseq.add(one, false), true);
        if decide(
            counts,
            |v| v == 0,
            name,
            "OIDLIST with MAXOID, LENGTH and NSEQ",
        )? {
            log.line(
                BRIEF,
                &[
                    prefix,
                    b"  [ERROR] OIDLIST cannot be provided without MAXOID , LENGTH, or NSEQ",
                ],
            )?;
            faults += 1;
        }
    }
    Ok(faults)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:794-804
/// The opening of the database of the alias file (`CSeqDB` of its path without the extension, as
/// a protein database for an extension `.p..`): 1 fault where NCBI's reader throws.
fn open(
    log: &mut Log<'_>,
    prefix: &[u8],
    file: &FileParts<'_>,
    path: &SearchPath,
) -> Result<u32, Stop> {
    let db = make_path(&file.dir, file.base);
    let seqtype = if file.ext.get(1) == Some(&b'p') {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    // The parts of a name that is a `String` cut at ASCII bytes: UTF-8 as the name.
    let db = String::from_utf8_lossy(&db);
    match SeqDb::open_in(&db, seqtype, path) {
        Ok(_) => Ok(0),
        Err(error) => {
            let exception = caught(&error)?;
            log.line(
                BRIEF,
                &[
                    prefix,
                    b"  [ERROR] caught exception in initializing blastdb",
                ],
            )?;
            log.what(prefix, exception.what.as_deref(), &|| unknown_text(&error))?;
            Ok(1)
        }
    }
}

#[cfg(test)]
#[path = "alias_tests.rs"]
mod tests;
