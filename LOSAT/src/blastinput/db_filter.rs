//! The taxid options of a database search (line DI-8): `-taxids`, `-taxidlist`,
//! `-negative_taxids`, `-negative_taxidlist` and `-no_taxid_expansion` of BLASTN, BLASTP,
//! BLASTX, TBLASTN and TBLASTX. The list is read and expanded as NCBI's handler of the
//! database arguments reads it (`s_GetTaxIDList`), and the search's database is restricted
//! to it (`seqdb::SeqDb::restrict_taxids`).
//!
//! `-seqidlist`, `-negative_seqidlist`, `-gilist` and `-negative_gilist` stay refused
//! (DI-8b); NCBI's exclusions with them are therefore not declared here.

use std::collections::BTreeSet;
use std::io::Write;

use anyhow::Result;
use clap::Args;

use crate::blastinput::app::options_error;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::seqdb::file::DbFile;
use crate::seqdb::tax_sqlite::TaxonomySqlite;
use crate::seqdb::{SearchPath, TaxidList};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2263-2299
// ```c++
//         arg_desc.AddOptionalKey(kArgTaxIdList, "taxids",
//                                 "Restrict search of database to include only "
//                                 "the specified taxonomy IDs and their descendants "
//                                 "(multiple IDs delimited by ',')",
//                                 CArgDescriptions::eString);
//         arg_desc.SetDependency(kArgTaxIdList, CArgDescriptions::eExcludes, kArgGiList);
//         arg_desc.SetDependency(kArgTaxIdList, CArgDescriptions::eExcludes, kArgSeqIdList);
//         arg_desc.SetDependency(kArgTaxIdList, CArgDescriptions::eExcludes, kArgNegativeGiList);
//         arg_desc.SetDependency(kArgTaxIdList, CArgDescriptions::eExcludes, kArgNegativeSeqidList);
//         arg_desc.SetDependency(kArgTaxIdList, CArgDescriptions::eExcludes, kArgRemote);
//         arg_desc.SetDependency(kArgTaxIdList, CArgDescriptions::eExcludes, kArgNegativeTaxIdList);
//         arg_desc.SetDependency(kArgTaxIdList, CArgDescriptions::eExcludes, kArgNegativeTaxIdListFile);
//
//         arg_desc.AddOptionalKey(kArgNegativeTaxIdList, "taxids",
//                                 "Restrict search of database to everything"
//                                 " except the specified taxonomy IDs and their descendants "
//                                 "(multiple IDs delimited by ',')",
//                                 CArgDescriptions::eString);
//         arg_desc.SetDependency(kArgNegativeTaxIdList, CArgDescriptions::eExcludes, kArgGiList);
//         arg_desc.SetDependency(kArgNegativeTaxIdList, CArgDescriptions::eExcludes, kArgSeqIdList);
//         arg_desc.SetDependency(kArgNegativeTaxIdList, CArgDescriptions::eExcludes, kArgNegativeGiList);
//         arg_desc.SetDependency(kArgNegativeTaxIdList, CArgDescriptions::eExcludes, kArgNegativeSeqidList);
//         arg_desc.SetDependency(kArgNegativeTaxIdList, CArgDescriptions::eExcludes, kArgRemote);
//
//         arg_desc.AddOptionalKey(kArgTaxIdListFile, "filename",
//                                 "Restrict search of database to include only "
//                                 "the specified taxonomy IDs and their descendants ",
//                                 CArgDescriptions::eString);
//         arg_desc.SetDependency(kArgTaxIdListFile, CArgDescriptions::eExcludes, kArgGiList);
//         arg_desc.SetDependency(kArgTaxIdListFile, CArgDescriptions::eExcludes, kArgSeqIdList);
//         arg_desc.SetDependency(kArgTaxIdListFile, CArgDescriptions::eExcludes, kArgNegativeGiList);
//         arg_desc.SetDependency(kArgTaxIdListFile, CArgDescriptions::eExcludes, kArgNegativeSeqidList);
//         arg_desc.SetDependency(kArgTaxIdListFile, CArgDescriptions::eExcludes, kArgRemote);
//         arg_desc.SetDependency(kArgTaxIdListFile, CArgDescriptions::eExcludes, kArgTaxIdList);
//         arg_desc.SetDependency(kArgTaxIdListFile, CArgDescriptions::eExcludes, kArgNegativeTaxIdList);
//         arg_desc.SetDependency(kArgTaxIdListFile, CArgDescriptions::eExcludes, kArgNegativeTaxIdListFile);
// ```
// The exclusions of the four lists with one another, with `-subject` and (on the
// programs' `-subject_loc`) with `-subject_loc` are clap's: LOSAT's parser reports them with
// its own text and exit status 2 (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, as `-db` with
// `-subject`; NCBI: usage, `Incompatible with argument`, exit 1).
/// The taxid options of a search program's command line.
#[derive(Args, Debug, Clone, Default)]
pub struct TaxidFilterArgs {
    /// `-taxids`: taxids separated by commas.
    #[arg(
        long = "taxids",
        value_name = "taxids",
        conflicts_with_all = ["taxidlist", "negative_taxids", "negative_taxidlist", "subject"],
        help = "Restrict search of database to include only the specified taxonomy IDs and their descendants (multiple IDs delimited by ',')"
    )]
    pub taxids: Option<String>,
    /// `-taxidlist`: a file of taxids, one per line.
    #[arg(
        long = "taxidlist",
        value_name = "filename",
        conflicts_with_all = ["negative_taxids", "negative_taxidlist", "subject"],
        help = "Restrict search of database to include only the specified taxonomy IDs and their descendants"
    )]
    pub taxidlist: Option<String>,
    /// `-negative_taxids`: taxids separated by commas.
    #[arg(
        long = "negative_taxids",
        value_name = "taxids",
        conflicts_with_all = ["negative_taxidlist", "subject"],
        help = "Restrict search of database to everything except the specified taxonomy IDs and their descendants (multiple IDs delimited by ',')"
    )]
    pub negative_taxids: Option<String>,
    /// `-negative_taxidlist`: a file of taxids, one per line.
    #[arg(
        long = "negative_taxidlist",
        value_name = "filename",
        conflicts_with_all = ["subject"],
        help = "Restrict search of database to everything except the specified taxonomy IDs and their descendants"
    )]
    pub negative_taxidlist: Option<String>,
    /// `-no_taxid_expansion`: the taxids alone, without their descendants.
    #[arg(
        long = "no_taxid_expansion",
        conflicts_with_all = ["subject"],
        help = "Do not expand the taxonomy IDs provided to their descendant taxonomy IDs "
    )]
    pub no_taxid_expansion: bool,
}

/// The names of the options of `TaxidFilterArgs`, for a program's `-subject_loc` exclusions.
pub const TAXID_OPTIONS: [&str; 5] = [
    "taxids",
    "taxidlist",
    "negative_taxids",
    "negative_taxidlist",
    "no_taxid_expansion",
];

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2428-2431
// ```c++
//         catch(CException &){
//             LOG_POST(Warning << "The -taxids command line option requires additional data files. Please see the section 'Taxonomic filtering for BLAST databases' in https://www.ncbi.nlm.nih.gov/books/NBK569839/ for details.");
//         }
// ```
const NO_TAXONOMY_WARNING: &str = "The -taxids command line option requires additional data \
     files. Please see the section 'Taxonomic filtering for BLAST databases' in \
     https://www.ncbi.nlm.nih.gov/books/NBK569839/ for details.\n";

/// NCBI's default name of the taxonomy file (`kDefaultName`, tax4blastsqlite.cpp).
const TAXONOMY_FILE: &str = "taxonomy4blast.sqlite3";

impl TaxidFilterArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2492-2503
    // ```c++
    //         } else if (args.Exist(kArgTaxIdList) && args[kArgTaxIdList]) {
    //         	s_GetTaxIDList(args[kArgTaxIdList].AsString(), false, false, m_SearchDb,args[kArgNoTaxIdExpansion].AsBoolean());
    //
    //         } else if (args.Exist(kArgTaxIdListFile) && args[kArgTaxIdListFile]) {
    //         	s_GetTaxIDList(args[kArgTaxIdListFile].AsString(), true, false, m_SearchDb, args[kArgNoTaxIdExpansion].AsBoolean());
    //
    //         } else if (args.Exist(kArgNegativeTaxIdList) && args[kArgNegativeTaxIdList]) {
    //         	s_GetTaxIDList(args[kArgNegativeTaxIdList].AsString(), false, true, m_SearchDb, args[kArgNoTaxIdExpansion].AsBoolean());
    //
    //         } else if (args.Exist(kArgNegativeTaxIdListFile) && args[kArgNegativeTaxIdListFile]) {
    //         	s_GetTaxIDList(args[kArgNegativeTaxIdListFile].AsString(), true, true, m_SearchDb,args[kArgNoTaxIdExpansion].AsBoolean());
    // ```
    /// Whether the command line gives one of the four lists.
    pub fn given(&self) -> bool {
        self.list().is_some()
    }

    /// The list option that NCBI reads: `(value, is a file, is negative)`.
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
    // ```c++
    // static void s_GetTaxIDList(const string & in, bool isFile, bool isNegativeList, CRef<CSearchDatabase> & sdb, bool isTargetOnly )
    // {
    //     vector<string> ids;
    //     if (isFile) {
    //         string filename(SeqDB_ResolveDbPath(in));
    //         if(filename == kEmptyStr) {
    //             NCBI_THROW(CInputException, eInvalidInput, "File is not acessible: "+ in );
    //         }
    //         CNcbiIfstream instream(filename.c_str());
    //         CStreamLineReader reader(instream);
    //
    //         while (!reader.AtEOF()) {
    //             reader.ReadLine();
    //             ids.push_back(reader.GetCurrentLine());
    //         }
    //     } else {
    //         NStr::Split(in, ",", ids, NStr::fSplit_Tokenize);
    //     }
    //     unique_ptr<ITaxonomy4Blast> tb;
    //     if( !isTargetOnly ) {
    //         try{
    //             tb.reset(new CTaxonomy4BlastSQLite());
    //         }
    //         catch(CException &){
    //             LOG_POST(Warning << "The -taxids command line option requires additional data files. Please see the section 'Taxonomic filtering for BLAST databases' in https://www.ncbi.nlm.nih.gov/books/NBK569839/ for details.");
    //         }
    //     }
    //     set<TTaxId> tax_ids;
    //     for (auto id : ids) {
    //         try {
    //             if (NStr::IsBlank(id)) {
    //                 continue;
    //             }
    //             auto taxid = NStr::StringToNumeric<TTaxId>(id, NStr::fAllowLeadingSpaces | NStr::fAllowTrailingSpaces);
    //             if( isTargetOnly ) {
    //                 tax_ids.insert(taxid);
    //             } else if (tb) {
    //                 tax_ids.insert(taxid);
    //                 vector<int> desc;
    //                 tb->GetLeafNodeTaxids(taxid, desc);
    //                 for (auto i: desc)
    //                     tax_ids.insert( static_cast<TTaxId>(i) );
    //             }
    //         } catch(CException &){
    //             NCBI_THROW(CInputException, eInvalidInput, "Invalid taxidlist file ");
    //         }
    //     }
    // ```
    // Without the taxonomy file NCBI reads every id and keeps none: the set is empty and the
    // search runs on the whole database (DI-8 DECISIONS 2 and 11, the Owner's decision of
    // 2026-10-29). A taxid whose descendants loop back to it (taxid 1, whose parent is 1)
    // never ends NCBI's recursive query; LOSAT visits each taxid once and ends with every
    // taxid below it (DECISIONS 10, `TaxonomySqlite::descendants`).
    /// The list of the command line, read and expanded as NCBI reads it, with the warning
    /// NCBI writes when the taxonomy file is missing written to `warnings`. `Ok(None)`: no
    /// list, or one without a taxid (NCBI then searches the whole database).
    pub fn resolve(
        &self,
        settings: &ApplicationSettings,
        warnings: &mut dyn Write,
    ) -> Result<Option<TaxidList>> {
        let Some((value, is_file, negative)) = self.list() else {
            return Ok(None);
        };
        // The registry's errors are the database's: `open_search_database` reports them.
        let registry = settings.blastdb.clone().unwrap_or(None);
        let path = SearchPath::from_environment(registry.as_deref())
            .map_err(|error| anyhow::anyhow!("{error}"))?;
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
                Err(_) => {
                    warnings.write_all(NO_TAXONOMY_WARNING.as_bytes())?;
                    None
                }
            },
            None => None,
        };
        let invalid = || options_error("Invalid taxidlist file ");
        let mut tax_ids = BTreeSet::new();
        for id in &ids {
            if id.iter().all(u8::is_ascii_whitespace) {
                continue;
            }
            let taxid = string_to_int(id).ok_or_else(invalid)?;
            if self.no_taxid_expansion {
                tax_ids.insert(taxid);
            } else if let Some(taxonomy) = &taxonomy {
                tax_ids.insert(taxid);
                tax_ids.extend(taxonomy.descendants(taxid).map_err(|_| invalid())?);
            }
        }
        if tax_ids.is_empty() {
            return Ok(None);
        }
        Ok(Some(if negative {
            TaxidList::Negative(tax_ids)
        } else {
            TaxidList::Positive(tax_ids)
        }))
    }
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:152-172
// ```c++
//     switch (m_EOLStyle) {
//     case eEOL_unknown: x_AdvanceEOLUnknown();                   break;
//     case eEOL_cr:      x_AdvanceEOLSimple('\r', '\n');          break;
//     case eEOL_lf:      x_AdvanceEOLSimple('\n', '\r');          break;
//     case eEOL_crlf:    x_AdvanceEOLCRLF();                      break;
//     case eEOL_mixed:   NcbiGetline(*m_Stream, m_Line, "\r\n");  break;
//     }
// ```
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

// `NStr::StringToNumeric<int>` with `fAllowLeadingSpaces | fAllowTrailingSpaces`
// (corelib/ncbistr.cpp, `NStr::StringToInt`): white space may surround the number; a sign (`+` or `-`) may lead it; the digits are
// decimal and the value must fit a signed 32-bit integer.
/// `NStr::StringToNumeric<TTaxId>` of a list entry, `None` for NCBI's exception.
fn string_to_int(token: &[u8]) -> Option<i32> {
    let text = std::str::from_utf8(token).ok()?;
    let text = text.trim_matches(|c: char| c.is_ascii_whitespace());
    let (negative, digits) = match text.as_bytes().first()? {
        b'-' => (true, &text[1..]),
        b'+' => (false, &text[1..]),
        _ => (false, text),
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let magnitude: i64 = digits.parse().ok()?;
    i32::try_from(if negative { -magnitude } else { magnitude }).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_end_as_ncbis_stream_reader_ends_them() {
        let lines = |text: &str| -> Vec<String> {
            stream_lines(text.as_bytes())
                .into_iter()
                .map(|line| String::from_utf8(line).unwrap())
                .collect()
        };
        assert_eq!(lines("9606\n10090\n"), ["9606", "10090"]);
        assert_eq!(
            lines("9606\r\n\r\n  10090  \r\n\n"),
            ["9606", "", "  10090  ", ""]
        );
        assert_eq!(lines("9606"), ["9606"]);
        assert_eq!(lines("1\r2"), ["1", "2"]);
        assert!(lines("").is_empty());
    }

    #[test]
    fn numbers_are_read_as_ncbis_string_to_numeric_reads_them() {
        assert_eq!(string_to_int(b"9606"), Some(9606));
        assert_eq!(string_to_int(b"+9606"), Some(9606));
        assert_eq!(string_to_int(b"  -5 "), Some(-5));
        assert_eq!(string_to_int(b"2147483647"), Some(i32::MAX));
        assert_eq!(string_to_int(b"-2147483648"), Some(i32::MIN));
        for bad in [
            &b"abc"[..],
            b"# human",
            b"2147483648",
            b"9606.0",
            b"+",
            b"1 2",
            b"0x10",
        ] {
            assert_eq!(string_to_int(bad), None, "{}", String::from_utf8_lossy(bad));
        }
    }
}
