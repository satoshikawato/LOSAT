//! The taxonomy of `losat blastdbcmd`: the taxids of `-taxids` / `-taxidlist` with their
//! descendants (`x_InitBlastDB_TaxIdList`), whose sequences `retrieve` writes, and
//! `-tax_info`, the taxids of the database (`x_PrintBlastDatabaseTaxInformation`).

use std::collections::BTreeSet;
use std::io::{BufWriter, Write};

use crate::blastinput::db_filter::{string_to_int, NO_TAXONOMY_WARNING, TAXONOMY_FILE};
use crate::seqdb::file::DbFile;
use crate::seqdb::tax_sqlite::TaxonomySqlite;
use crate::seqdb::taxdb::{LazyTaxDb, TaxNames};
use crate::seqdb::{SearchPath, SeqDb};

use super::args::{self, Args};
use super::{database_error, refusal, Diag};

/// What ends the reading of the taxids before the database opens: the exit status, the
/// diagnostics already posted.
pub(super) struct Ended(pub i32);

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:320-357
// ```c++
// void
// CBlastDBCmdApp::x_InitBlastDB_TaxIdList()
// {
//    	const CArgs& args = GetArgs();
//     vector<string> ids;
//     if(args[kArgTaxIdList].HasValue()) {
//     	string input = args[kArgTaxIdList].AsString();
//     	NStr::Split(input, ",", ids);
//     }
//     else {
//     	CNcbiIstream& input = args[kArgTaxIdListFile].AsInputFile();
//     	while (input) {
//     		string line;
//     	    NcbiGetlineEOL(input, line);
//     	    if ( !line.empty() ) {
//     	       	ids.push_back(line);
//     	    }
//     	}
//     }
//
//     unique_ptr<ITaxonomy4Blast> tb;
//     if( ! args[kArgNoTaxIdExpansion].AsBoolean() ){
//         try{
//             tb.reset(new CTaxonomy4BlastSQLite());
//         }
//         catch(CException &){
//             LOG_POST(Warning << "The -taxids command line option requires additional data files. Please see the section 'Taxonomic filtering for BLAST databases' in https://www.ncbi.nlm.nih.gov/books/NBK569839/ for details.");
//         }
//     }
//
//     for (auto id : ids) {
//         auto taxid = NStr::StringToNumeric<TTaxId>(id, NStr::fAllowLeadingSpaces | NStr::fAllowTrailingSpaces);
//     	m_TaxIdList.insert(taxid);
//         if (tb) {
//             vector<int> descendants;
//             tb->GetLeafNodeTaxids(taxid, descendants);
//             for (auto d: descendants)
//                 m_TaxIdList.insert(static_cast<TTaxId>(d));
//         }
//     }
// ```
// `NStr::Split` keeps the empty pieces of `1,,2` and `1,` and gives none for the empty
// string (probes J101, J102, S36 of line BDC). The lines of the file are those of
// `NcbiGetlineEOL`, which on Linux ends a line at `\n` alone (a `\r` stays and reads as white
// space around the number: J109), as `-entry_batch` reads them. Unlike a search program's
// list (`db_filter`), the id is kept without the taxonomy file (T02), and a line of blanks
// is read as a number. An id that is not a number escapes NCBI's `Run` as a
// `CStringException` with the path of NCBI's build in its text (exit status 255: J101,
// J102, J107, J108, J131, S35), which LOSAT refuses (line BDC DECISIONS 8). `-taxids 1`
// never ends in NCBI (the parent of taxid 1 is 1, probe J120); LOSAT's walk visits each
// taxid once and ends with every taxid below it (`TaxonomySqlite::descendants`, DI-8
// DECISIONS 10).
/// The taxids of `-taxids` / `-taxidlist`, with their descendants in the taxonomy file of
/// the search path unless `-no_taxid_expansion`; NCBI's warning when that file is missing
/// or not one NCBI reads goes to `diag`.
pub(super) fn tax_id_list(
    args: &Args,
    path: &SearchPath,
    diag: &mut Diag,
) -> Result<BTreeSet<i32>, Ended> {
    let ids: Vec<Vec<u8>> = match (args.value("taxids"), args.value("taxidlist")) {
        (Some(""), _) => Vec::new(),
        (Some(input), _) => input.split(',').map(|id| id.as_bytes().to_vec()).collect(),
        (None, Some(file)) => {
            let Ok(bytes) = args::read_input_file("taxidlist", file) else {
                // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
                // ```c++
                //     catch (const CArgException& e) {                                        \
                //         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
                //         exit_code = BLAST_INPUT_ERROR;                                      \
                //     }                                                                       \
                // ```
                diag.post(&format!(
                    "Command line argument error: {}",
                    args::arg_message("taxidlist", "File is not accessible", file)
                ));
                return Err(Ended(1));
            };
            bytes
                .split(|&byte| byte == b'\n')
                .filter(|line| !line.is_empty())
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
            Err(_) => {
                diag.post(NO_TAXONOMY_WARNING.trim_end_matches('\n'));
                None
            }
        },
        None => None,
    };
    let mut tax_ids = BTreeSet::new();
    for id in &ids {
        let Some(taxid) = string_to_int(id) else {
            return Err(Ended(diag.refuse(&refusal(&format!(
                "the taxid '{}' is not a number: NCBI ends with an exception that names its \
                 source files (exit status 255)",
                String::from_utf8_lossy(id)
            )))));
        };
        tax_ids.insert(taxid);
        if let Some(taxonomy) = &taxonomy {
            let descendants = taxonomy.descendants(taxid).map_err(|reason| {
                Ended(diag.refuse(&refusal(&format!(
                    "{TAXONOMY_FILE}: the descendants of taxid {taxid} are not readable \
                     ({reason})"
                ))))
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

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:592-647
// ```c++
// 	CPrintTaxFields(CNcbiOstream & out, const string & fmt): m_Out(out), m_NeedTaxInfoLookup(false), m_NeedNumSeqs(false) {
// 		vector<string> fields;
//         string sp = kEmptyStr;
// 	    for (unsigned int i = 0; i < fmt.size(); i++) {
// 	    	if (fmt[i] == '%') {
// 		        if (fmt[i+1] == '%') {
// 		            sp += fmt[i];
// 		            continue;
// 		        }
// 		        i++;
// 		        switch (fmt[i]) {
// 		        case 'f' :
// 		        	m_NeedTaxInfoLookup = true;
// 		      		m_Seperators.push_back(sp);
// 		      		for(unsigned int i=eTaxID; i < eMaxFields; i++){
// 		                if ( i == eNumSeqs ) continue;
// 		      		    m_Fields.push_back(i);
// 		      		    m_Seperators.push_back("\t");
// 		      		}
// 		        break;
// 		        case 'T' :
//        				m_Fields.push_back(eTaxID);
//        			break;
// 		        ...
// 		        case 'n' :
// 		        	m_Fields.push_back(eNumSeqs);
//                     m_NeedNumSeqs = true;
// 		        break;
// 		        default:
// 	                sp += fmt[i-1];
// 	                sp += fmt[i];
// 	                continue;
// 	            break;
// 		        }
// 	            m_Seperators.push_back(sp);
// 	            sp = kEmptyStr;
// 	    	}
// 	    	else {
// 		        sp += fmt[i];
// 	    	}
// 	    }
// 		m_Seperators.push_back(sp);
//
// 		if(m_Fields.empty()) {
// 			NCBI_THROW(CInputException, eInvalidInput,
// 				       "Invalid format options for tax_info.");
// 		}
// 	}
// ```
// Ported as it reads, the quirks with it (probes T10-T12, J140-J159 of line BDC): `%f`
// pushes the text before it and a tab per field, then the text after the `switch` again,
// so after `%f` the fields shift onto the tabs (`%f %n` is `T\tS\tL\tK\tB\tN`); `%%` adds one
// `%` and the second `%` starts a specifier (`%%T` is `%` and the taxid, `%%x` the text
// `%%x`); a `%` that ends the format reads the string's NUL and adds `%` and a NUL byte;
// any other letter adds itself and its `%` as text.
impl TaxFields {
    /// `CPrintTaxFields(out, fmt)`; `Err` for NCBI's `Invalid format options for tax_info.`
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

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:665-697
    // ```c++
    // 	void PrintEntry(const SSeqDBTaxInfo & t, int num_seqs){
    // 		for(unsigned int i=0; i < m_Fields.size(); i++) {
    // 			m_Out << m_Seperators[i];
    // 			switch (m_Fields[i]){
    // 				case eTaxID:
    // 					m_Out << t.taxid;
    // 				break;
    // 				...
    // 				case eNumSeqs:
    // 					m_Out << num_seqs;
    // 				break;
    // 				...
    // 			}
    // 		}
    // 		m_Out << m_Seperators.back();
    // 		m_Out << "\n";
    // 	}
    // ```
    /// `PrintEntry`: the line of `taxid` with its `names` (`None`: NCBI's `N/A`) and its
    /// number of sequences.
    fn line(&self, taxid: i32, names: Option<&TaxNames>, num_seqs: usize, out: &mut Vec<u8>) {
        const NA: &[u8] = b"N/A";
        for (i, field) in self.fields.iter().enumerate() {
            out.extend_from_slice(&self.separators[i]);
            match field {
                Field::TaxId => out.extend_from_slice(taxid.to_string().as_bytes()),
                Field::SciName => out.extend_from_slice(names.map_or(NA, |n| &n.scientific_name)),
                Field::CommonName => out.extend_from_slice(names.map_or(NA, |n| &n.common_name)),
                Field::SuperKingdom => {
                    out.extend_from_slice(names.map_or(NA, |n| &n.super_kingdom))
                }
                Field::BlastName => out.extend_from_slice(names.map_or(NA, |n| &n.blast_name)),
                Field::NumSeqs => out.extend_from_slice(num_seqs.to_string().as_bytes()),
            }
        }
        out.extend_from_slice(self.separators.last().map_or(&[][..], Vec::as_slice));
        out.push(b'\n');
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:700-733
// ```c++
//     CNcbiOstream& out = args[kArgOutput].AsOutputFile();
//     const string& kFmt = args["outfmt"].AsString();
//     CPrintTaxFields tf(out, kFmt);
//     set<TTaxId> tax_ids;
//     m_BlastDb->GetDBTaxIds(tax_ids);
//
//     // Print basic database information
//     out << "# of Tax IDs in Database: " << tax_ids.size() << endl;
// 	SSeqDBTaxInfo info;
//     ITERATE(set<TTaxId>, itr, tax_ids) {
//     	SSeqDBTaxInfo info;
//     	if(tf.NeedTaxNames()){
//     		CSeqDBTaxInfo::GetTaxNames(*itr, info);
//     		if(info.taxid == ZERO_TAX_ID){
//     			info.taxid = *itr;
//    				info.scientific_name = NA;
//    				...
//    			}
//     	}
//     	else {
//    			info.taxid = *itr;
//     	}
//         int num_seqs = 0;
//         if(tf.NeedNumSeqs()) {
//             vector<blastdb::TOid> rv;
//             set<TTaxId> ti;
//             ti.insert(*itr);
//             m_BlastDb->TaxIdsToOids(ti, rv);
//             num_seqs = rv.size();
//         }
//    		tf.PrintEntry(info, num_seqs);
//     }
// ```
// NCBI's exceptions escape to `Run`: the format's (`CInputException`, exit status 1), the
// reader's (`CSeqDBException`, exit status 2) after the lines already written. Taxid 0 (the
// sequences without one) is a taxid like any other; taxdb has no names for it (T12).
/// `-tax_info` of `db`: the number of its taxids, then a line for each, ascending, to `out`;
/// the exit status.
pub(super) fn tax_info(
    db: &SeqDb,
    args: &Args,
    path: &SearchPath,
    diag: &mut Diag,
    out: &mut dyn Write,
) -> i32 {
    let fmt = args.value("outfmt").unwrap_or("%f");
    let Ok(fields) = TaxFields::parse(fmt.as_bytes()) else {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-176
        // ```c++
        //     catch (const blast::CInputException& e) {                               \
        //         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
        //         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
        //         exit_code = BLAST_INPUT_ERROR;                                      \
        //     }                                                                       \
        // ```
        diag.post("BLAST query/options error: Invalid format options for tax_info.");
        diag.post("Please refer to the BLAST+ user manual.");
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
    let taxdb = LazyTaxDb::new(path.clone());
    let mut sink = BufWriter::new(out);
    let _ = writeln!(sink, "# of Tax IDs in Database: {}", tax_ids.len());
    let mut line = Vec::new();
    for &taxid in &tax_ids {
        let names = if fields.need_names {
            let mut messages = Vec::new();
            match taxdb.names(taxid, &mut messages) {
                Ok(names) if messages.is_empty() => names,
                Ok(_) => {
                    let _ = sink.flush();
                    return diag.refuse(&refusal(&format!(
                        "a taxonomy database NCBI opens with a diagnostic ({})",
                        messages[0]
                    )));
                }
                Err(error) => {
                    let _ = sink.flush();
                    return database_error(diag, error);
                }
            }
        } else {
            None
        };
        let num_seqs = if fields.need_num_seqs {
            match lookup.tax_ids_to_oids(&BTreeSet::from([taxid])) {
                Ok(oids) => oids.len(),
                Err(error) => {
                    let _ = sink.flush();
                    return database_error(diag, error);
                }
            }
        } else {
            0
        };
        line.clear();
        fields.line(taxid, names.as_ref(), num_seqs, &mut line);
        // NCBI's stream ignores a failed write (no exception is enabled on it).
        let _ = sink.write_all(&line);
    }
    let _ = sink.flush();
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(fmt: &str) -> Vec<u8> {
        let fields = TaxFields::parse(fmt.as_bytes()).unwrap();
        let names = TaxNames {
            taxid: 9606,
            scientific_name: b"Homo sapiens".to_vec(),
            common_name: b"human".to_vec(),
            blast_name: b"primates".to_vec(),
            super_kingdom: b"Eukaryota".to_vec(),
        };
        let mut out = Vec::new();
        fields.line(9606, Some(&names), 3, &mut out);
        out
    }

    // NCBI's lines for these formats: probes T10-T12, J140-J159 of line BDC.
    #[test]
    fn the_format_is_read_as_ncbi_s_loop_reads_it() {
        assert_eq!(
            line("%f"),
            b"9606\tHomo sapiens\thuman\tEukaryota\tprimates\n"
        );
        assert_eq!(
            line("%f %n"),
            b"9606\tHomo sapiens\thuman\tEukaryota\tprimates\t3\n"
        );
        assert_eq!(
            line("%T %f"),
            b"9606 9606\tHomo sapiens\thuman\tEukaryota\tprimates\n"
        );
        assert_eq!(
            line("%T|%S|%L|%K|%B|%n"),
            b"9606|Homo sapiens|human|Eukaryota|primates|3\n"
        );
        assert_eq!(line("%T%n%S"), b"96063Homo sapiens\n");
        assert_eq!(line("%%T"), b"%9606\n");
        assert_eq!(line("%%x %T"), b"%%x 9606\n");
        assert_eq!(line("%x %T %n"), b"%x 9606 3\n");
        assert_eq!(line("%T %"), b"9606 %\0\n");
        assert_eq!(line("%T %%"), b"9606 %%\0\n");
        assert_eq!(line("%n"), b"3\n");
        assert!(TaxFields::parse(b"abc").is_err());
        assert!(TaxFields::parse(b"").is_err());
        assert!(TaxFields::parse(b"%x").is_err());
    }

    #[test]
    fn missing_names_are_not_available() {
        let fields = TaxFields::parse(b"%f").unwrap();
        let mut out = Vec::new();
        fields.line(0, None, 0, &mut out);
        assert_eq!(out, b"0\tN/A\tN/A\tN/A\tN/A\n");
    }
}
