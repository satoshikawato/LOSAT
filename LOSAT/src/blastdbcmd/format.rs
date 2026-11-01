//! `-outfmt` with `%` specifiers other than `%f`: NCBI's `CBlastDB_SeqFormatter` (one line
//! per definition line, or one per sequence when no specifier reads a definition line) and
//! the fields `CBlastDeflineUtil::ExtractDataFromBlastDefline` extracts.

use std::collections::BTreeSet;
use std::io::Write;

use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::label::{bare, fasta};
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::{BlastDefline, DeflineSet, ParsedSeqId};

use super::letters::{hash_text, sequence_hash};
use super::retrieve::{numeric_target, parse_target, Config, Source, Stop};

/// `NOT_AVAILABLE`.
const NOT_AVAILABLE: &str = "N/A";
/// `SEPARATOR`.
const SEPARATOR: &str = ";";

/// The `-outfmt` text read into its specifiers (`m_ReplTypes`) and the literal text
/// between them (`m_Seperators`, one more than the specifiers).
#[derive(Debug)]
pub(super) struct FormatSpec {
    types: Vec<u8>,
    separators: Vec<Vec<u8>>,
    /// `m_GetDefline`: a specifier reads the definition lines.
    get_defline: bool,
    /// `m_OtherFields & (1 << e_seq)`, `(1 << e_hash)`.
    need_seq: bool,
    need_hash: bool,
}

impl FormatSpec {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:59-91
    // ```c++
    //     // Record where the offsets where the replacements must occur
    // 	string sp = kEmptyStr;
    //     for (SIZE_TYPE i = 0; i < m_FmtSpec.size(); i++) {
    //         if (m_FmtSpec[i] == '%') {
    //         	if ( m_FmtSpec[i+1] == '%') {
    //         	    // remove the escape character for '%'
    //         		i++;
    //         		sp += m_FmtSpec[i];
    //         	    continue;
    //         	}
    //         	i++;
    //             m_ReplTypes.push_back(m_FmtSpec[i]);
    //             m_Seperators.push_back(sp);
    //             sp = kEmptyStr;
    //         }
    //         else {
    //         	sp += m_FmtSpec[i];
    //         }
    //     }
    //     m_Seperators.push_back(sp);
    //
    //     if (m_ReplTypes.empty() || (m_ReplTypes.size() + 1 != m_Seperators.size())) {
    //         NCBI_THROW(CInvalidDataException, eInvalidInput,
    //                    "Invalid format specification");
    //     }
    // ```
    /// Reads the format `spec` (the `-outfmt` value after blastdbcmd's handling of `%f`,
    /// `%b`, `%d` and `%m`); `Err` is the message of NCBI's exception. A `%` that ends the
    /// text reads the string's terminating NUL as its specifier, as NCBI's loop does.
    pub fn parse(spec: &str) -> Result<Self, String> {
        let bytes = spec.as_bytes();
        let mut types = Vec::new();
        let mut separators = Vec::new();
        let mut sp = Vec::new();
        let mut i = 0;
        while i < bytes.len() {
            if bytes[i] == b'%' {
                if bytes.get(i + 1) == Some(&b'%') {
                    i += 1;
                    sp.push(b'%');
                } else {
                    i += 1;
                    types.push(bytes.get(i).copied().unwrap_or(0));
                    separators.push(std::mem::take(&mut sp));
                }
            } else {
                sp.push(bytes[i]);
            }
            i += 1;
        }
        separators.push(sp);
        if types.is_empty() {
            return Err("Invalid format specification".to_string());
        }
        let get_defline = types.iter().any(|t| b"aigtTXPLBKSCNend".contains(t));
        Ok(FormatSpec {
            need_seq: types.contains(&b's'),
            need_hash: types.contains(&b'h'),
            types,
            separators,
            get_defline,
        })
    }

    /// Whether a specifier writes the letters of the sequence (`%s`).
    pub fn writes_letters(&self) -> bool {
        self.need_seq
    }

    /// Whether a specifier asks for the ASN.1 text of a definition line (`%d`), which LOSAT
    /// does not write (line BDC DECISIONS 9).
    pub fn asks_asn_text(&self) -> bool {
        self.types.contains(&b'd')
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:315-358
    // ```c++
    // int CBlastDB_SeqFormatter::Write(CSeqDB::TOID oid, const CBlastDB_FormatterConfig & config, string target_id)
    // {
    // 	int status = -1;
    //     vector<string> other_fields(e_max_other_fields, kEmptyStr);
    //     string mask = kEmptyStr;
    //     if (m_OtherFields) {
    //     	if (m_OtherFields & (1 << e_seq)) {
    //     		x_GetSeq(oid, config, other_fields[e_seq]);
    //     	}
    //     	if (m_OtherFields & (1 << e_hash)) {
    //     		other_fields[e_hash] = x_GetSeqHash(oid);
    //     	}
    //     	...
    //     }
    //     if (m_GetDefline) {
    //     	CRef<CBlast_def_line_set> df_set = m_BlastDb.GetHdr(oid);
    //     	if(df_set.Empty()) {
    //     		return status;
    //     	}
    //     	if (target_id == kEmptyStr) {
    //     		...
    //     		ITERATE(CBlast_def_line_set::Tdata, itr, df_set->Get()) {
    //     			vector<string> defline_data(CBlastDeflineUtil::max_index, kEmptyStr);
    //     			CBlastDeflineUtil::ExtractDataFromBlastDefline(**itr, defline_data, m_DeflineFields, false);
    //     			x_Print(oid, defline_data, other_fields);
    //     		}
    //     	}
    //     	else {
    //    			vector<string> defline_data(CBlastDeflineUtil::max_index, kEmptyStr);
    //     		CBlastDeflineUtil::ExtractDataFromBlastDeflineSet(*df_set, defline_data, m_DeflineFields, target_id, false);
    //    			x_Print(oid, defline_data, other_fields);
    //     	}
    //     }
    //     else {
    // 		vector<string> defline_data(CBlastDeflineUtil::max_index, kEmptyStr);
    // 		x_Print(oid, defline_data, other_fields);
    //     }
    //     return 0;
    // }
    // ```
    /// Writes the lines of sequence `oid` to `out`: one per definition line (the one of
    /// `target` alone with `-target_only`), or one line when no specifier reads them. The
    /// line NCBI has written up to an unknown specifier is written before its error.
    pub fn write(
        &self,
        source: &Source,
        oid: u32,
        config: &Config,
        target: Option<&str>,
        out: &mut dyn Write,
    ) -> Result<(), Stop> {
        let mut seq = Vec::new();
        let mut hash = String::new();
        if self.need_seq {
            seq = x_get_seq(source, oid, config)?;
        }
        if self.need_hash {
            hash = hash_text(sequence_hash(&source.string_letters(oid, None)?));
        }
        let length = source.length(oid)?;
        let mut line = Vec::new();
        let mut print = |defline: Option<&BlastDefline>, line: &mut Vec<u8>| -> Result<(), Stop> {
            line.clear();
            let printed = self.x_print(source, oid, length, defline, &seq, &hash, line);
            out.write_all(line).map_err(Stop::io)?;
            printed
        };
        if !self.get_defline {
            return print(None, &mut line);
        }
        let set = source.deflines(oid)?;
        if set.deflines.is_empty() {
            return Ok(());
        }
        match target {
            None => {
                for defline in &set.deflines {
                    print(Some(defline), &mut line)?;
                }
                Ok(())
            }
            Some(target) => {
                let defline = extract_target_defline(&set, target)?;
                print(Some(defline), &mut line)
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:187-201
    // ```c++
    //         case 'o':
    //             m_Out << NStr::NumericToString(oid);
    //             break;
    //         ...
    //         case 'l':
    //             m_Out << NStr::NumericToString(m_BlastDb.GetSeqLength(oid));
    //             break;
    // ```
    // `Write(-1, config)` (`-ipg` of a PIG no volume holds): the sequence, its hash and the
    // definition lines are read before the line (`x_GetSeq`, `x_GetSeqHash`, `GetHdr`), and
    // each read throws `OID not found`; without them the line is written up to `%l`, whose
    // length throws, and `%o` is `-1`.
    /// NCBI's `Write` of the OID -1: the part of the line it writes, then `OID not found`
    /// where it reads the database.
    pub fn write_missing_oid(&self, out: &mut dyn Write) -> Result<(), Stop> {
        let missing = || Stop::Ncbi(super::retrieve::OID_NOT_FOUND.to_string());
        if self.need_seq || self.need_hash || self.get_defline {
            return Err(missing());
        }
        let mut line = Vec::new();
        for (i, &kind) in self.types.iter().enumerate() {
            line.extend_from_slice(&self.separators[i]);
            match kind {
                b'o' => line.extend_from_slice(b"-1"),
                b'l' => {
                    out.write_all(&line).map_err(Stop::io)?;
                    return Err(missing());
                }
                b'm' | b'd' => {
                    return Err(Stop::Refused(format!(
                        "the -outfmt specifier %{} is not supported by LOSAT's {}",
                        kind as char,
                        super::PROGRAM
                    )));
                }
                other => {
                    out.write_all(&line).map_err(Stop::io)?;
                    let mut message = b"Unrecognized format specification: '%".to_vec();
                    message.push(other);
                    message.push(b'\'');
                    return Err(Stop::NcbiBytes(message));
                }
            }
        }
        line.extend_from_slice(&self.separators[self.types.len()]);
        line.push(b'\n');
        out.write_all(&line).map_err(Stop::io)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:165-262
    // ```c++
    // void CBlastDB_SeqFormatter::x_Print(CSeqDB::TOID oid, vector<string> & defline_data, vector<string> & other_fields)
    // {
    //     for(unsigned int i =0; i < m_ReplTypes.size(); i++) {
    //     	m_Out << m_Seperators[i];
    //         switch (m_ReplTypes[i]) {
    //
    //         case 's':
    //             m_Out << other_fields[e_seq];
    //             break;
    //         ...
    //         default:
    //             CNcbiOstrstream os;
    //             os << "Unrecognized format specification: '%" << m_ReplTypes[i] << "'";
    //             NCBI_THROW(CInvalidDataException, eInvalidInput,
    //                        CNcbiOstrstreamToString(os));
    //         }
    //     }
    //     m_Out << m_Seperators.back();
    //     m_Out << endl;
    // }
    // ```
    /// One line into `line`; `Err` after the part NCBI writes before an unknown specifier.
    fn x_print(
        &self,
        source: &Source,
        oid: u32,
        length: u32,
        defline: Option<&BlastDefline>,
        seq: &[u8],
        hash: &str,
        line: &mut Vec<u8>,
    ) -> Result<(), Stop> {
        let fields = match defline {
            Some(defline) => Some(DeflineFields::extract(source, defline, &self.types)?),
            None => None,
        };
        let empty = DeflineFields::default();
        let fields = fields.as_ref().unwrap_or(&empty);
        for (i, &kind) in self.types.iter().enumerate() {
            line.extend_from_slice(&self.separators[i]);
            match kind {
                b's' => line.extend_from_slice(seq),
                b'a' => line.extend_from_slice(&fields.accession),
                b'i' => line.extend_from_slice(&fields.seq_id),
                b'g' => line.extend_from_slice(fields.gi.as_bytes()),
                b'o' => line.extend_from_slice(source.db.oid(oid).to_string().as_bytes()),
                b't' => line.extend_from_slice(&fields.title),
                b'h' => line.extend_from_slice(hash.as_bytes()),
                b'l' => line.extend_from_slice(length.to_string().as_bytes()),
                b'T' => line.extend_from_slice(fields.tax_id.as_bytes()),
                b'X' => line.extend_from_slice(fields.leaf_node_tax_ids.as_bytes()),
                b'P' => line.extend_from_slice(fields.pig.as_bytes()),
                b'L' => line.extend_from_slice(&fields.common_name),
                b'C' => line.extend_from_slice(&fields.leaf_node_common_names),
                b'B' => line.extend_from_slice(&fields.blast_name),
                b'K' => line.extend_from_slice(&fields.super_kingdom),
                b'S' => line.extend_from_slice(&fields.scientific_name),
                b'N' => line.extend_from_slice(&fields.leaf_node_scientific_names),
                b'e' => line.extend_from_slice(fields.membership.as_bytes()),
                b'n' => line.extend_from_slice(fields.links.as_bytes()),
                b'm' | b'd' => {
                    // `%m` never reaches the formatter (blastdbcmd refuses the mask
                    // algorithm of a database without masks first); `%d` is refused before
                    // the first line.
                    return Err(Stop::Refused(format!(
                        "the -outfmt specifier %{} is not supported by LOSAT's {}",
                        kind as char,
                        super::PROGRAM
                    )));
                }
                other => {
                    let mut message = b"Unrecognized format specification: '%".to_vec();
                    message.push(other);
                    message.push(b'\'');
                    return Err(Stop::NcbiBytes(message));
                }
            }
        }
        line.extend_from_slice(&self.separators[self.types.len()]);
        line.push(b'\n');
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:282-313
// ```c++
// void CBlastDB_SeqFormatter::x_GetSeq(CSeqDB::TOID oid, const CBlastDB_FormatterConfig & config, string & seq)
// {
// 	TSeqRange r;
// 	if(config.m_SeqRange.NotEmpty()) {
// 		TSeqPos length = m_BlastDb.GetSeqLength(oid);
// 	    r = config.m_SeqRange;
// 	    if((TSeqPos)length <= config.m_SeqRange.GetTo()) {
// 	    	r.SetTo(length-1);
// 	    }
// 	}
// 	if(r.Empty()) {
// 	   	m_BlastDb.GetSequenceAsString(oid, seq);
// 	}
// 	else {
// 	   	m_BlastDb.GetSequenceAsString(oid, seq, r);
// 	}
// 	...
// 	if(config.m_Strand == eNa_strand_minus) {
// 		CBlastSeqUtil::GetReverseStrandSeq(seq);
// 	}
// }
// ```
/// The letters of `%s`: the range's (the whole sequence when the range is empty once
/// clipped), reverse-complemented for `-strand minus` (nucleotide databases: `Config` never
/// holds the minus strand for a protein one).
fn x_get_seq(source: &Source, oid: u32, config: &Config) -> Result<Vec<u8>, Stop> {
    let range = config.range.clipped(source.length(oid)?);
    let mut seq = source.string_letters(oid, range)?;
    if config.strand_minus {
        super::letters::reverse_complement(&mut seq);
    }
    Ok(seq)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:800-819
// ```c++
// void CBlastDeflineUtil::ExtractDataFromBlastDeflineSet(const CBlast_def_line_set & dl_set,
//         											   vector<string> & results,
//         											   BlastDeflineFields fields,
//         											   string target_id,
//         											   bool use_long_id)
// {
// 	CSeq_id target_seq_id (target_id, CSeq_id::fParse_PartialOK | CSeq_id::fParse_Default);
// 	Int8 num_id = NStr::StringToNumeric<Int8>(target_id, NStr::fConvErr_NoThrow);
// 	bool can_be_gi = errno ? false: true;
// 	ITERATE(CBlast_def_line_set::Tdata, itr, dl_set.Get()) {
// 		 ITERATE(CBlast_def_line::TSeqid, id, (*itr)->GetSeqid()) {
// 			 if ((*id)->Match(target_seq_id) || (can_be_gi && (*id)->IsGi() && ((*id)->GetGi() == GI_FROM(TIntId, num_id)))) {
// 				 CBlastDeflineUtil::ExtractDataFromBlastDefline( **itr, results, fields, use_long_id);
// 				 return;
// 			 }
// 		 }
// 	}
//
// 	NCBI_THROW(CException, eInvalid, "Failed to find target id " + target_id);
// }
// ```
/// The definition line of `-target_only` for the `%` formats: the first whose Seq-id
/// matches the target (`CSeq_id::Match`) or whose GI is the target's number.
fn extract_target_defline<'a>(set: &'a DeflineSet, target: &str) -> Result<&'a BlastDefline, Stop> {
    let target_id = parse_target(target)?;
    let number = numeric_target(target);
    for defline in &set.deflines {
        for id in &defline.ids {
            let matches = compare(id, &target_id) == Sic::Yes
                || number.is_some_and(|number| id.gi() == Some(number));
            if matches {
                return Ok(defline);
            }
        }
    }
    Err(Stop::Ncbi(format!("Failed to find target id {target}")))
}

/// The strings of `CBlastDeflineUtil::ExtractDataFromBlastDefline` (`results`), those the
/// format asks for; the others empty.
#[derive(Default)]
struct DeflineFields {
    gi: String,
    seq_id: Vec<u8>,
    accession: Vec<u8>,
    title: Vec<u8>,
    tax_id: String,
    scientific_name: Vec<u8>,
    common_name: Vec<u8>,
    blast_name: Vec<u8>,
    super_kingdom: Vec<u8>,
    leaf_node_tax_ids: String,
    leaf_node_scientific_names: Vec<u8>,
    leaf_node_common_names: Vec<u8>,
    membership: String,
    pig: String,
    links: String,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:821-827
// ```c++
// static string s_CheckName(const string & name)
// {
//         if(name == "-") return NOT_AVAILABLE;
//         if(name == "unclassified") return NOT_AVAILABLE;
//
//         return name;
// }
// ```
fn check_name(name: &[u8]) -> Vec<u8> {
    if name == b"-" || name == b"unclassified" {
        NOT_AVAILABLE.as_bytes().to_vec()
    } else {
        name.to_vec()
    }
}

impl DeflineFields {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:829-957
    // ```c++
    // 	if (fields.gi == 1) {
    //          results[CBlastDeflineUtil::gi] = NOT_AVAILABLE;
    // 		 ITERATE(CBlast_def_line::TSeqid, id, dl.GetSeqid()) {
    // 			 if ((*id)->IsGi()) {
    // 				 TGi gi = (*id)->GetGi();
    // 				 results[CBlastDeflineUtil::gi] = NStr::NumericToString(gi);
    // 				 break;
    // 			 }
    // 	     }
    // 	}
    // 	if ((fields.accession == 1) || (fields.seq_id == 1)) {
    // 		 CRef<CSeq_id> theId = FindBestChoice(dl.GetSeqid(), CSeq_id::WorstRank);
    // 		 if(fields.seq_id == 1) {
    //         	 results[CBlastDeflineUtil::seq_id] = theId->AsFastaString();
    // 		 }
    // 		 if(fields.accession == 1) {
    // 			 results[CBlastDeflineUtil::accession] = GetBareId(*theId);
    // 		 }
    // 	}
    // 	if(fields.title == 1) {
    // 		if(dl.IsSetTitle()) {
    //         	 results[CBlastDeflineUtil::title] = dl.GetTitle();
    // 		}
    // 		else {
    //         	 results[CBlastDeflineUtil::title] = NOT_AVAILABLE;
    // 		}
    // 	}
    // 	if ((fields.tax_id == 1) || (fields.tax_names == 1)) {
    //         TTaxId tax_id = ZERO_TAX_ID;
    // 		if (dl.IsSetTaxid()) {
    // 			tax_id = dl.GetTaxid();
    // 		}
    //
    // 		if (fields.tax_id == 1) {
    // 			results[CBlastDeflineUtil::tax_id] = NStr::NumericToString(tax_id);
    // 		}
    //
    // 		if (fields.tax_names == 1) {
    // 			try {
    // 				SSeqDBTaxInfo taxinfo;
    // 			    CSeqDB::GetTaxInfo(tax_id, taxinfo);
    // 			    results[CBlastDeflineUtil::scientific_name] = taxinfo.scientific_name;
    // 			    results[CBlastDeflineUtil::common_name] = taxinfo.common_name;
    // 			    results[CBlastDeflineUtil::blast_name] = s_CheckName(taxinfo.blast_name);
    // 			    results[CBlastDeflineUtil::super_kingdom] = s_CheckName(taxinfo.s_kingdom);
    // 	        } catch (const CException&) {
    // 			    results[CBlastDeflineUtil::scientific_name] = NOT_AVAILABLE;
    // 			    results[CBlastDeflineUtil::common_name] = NOT_AVAILABLE;
    // 			    results[CBlastDeflineUtil::blast_name] = NOT_AVAILABLE;
    // 			    results[CBlastDeflineUtil::super_kingdom] = NOT_AVAILABLE;
    // 			}
    // 		}
    // 	}
    //
    // 	if ((fields.leaf_node_tax_ids == 1) || (fields.leaf_node_tax_names == 1)) {
    // 		set<TTaxId>  tax_id_set = dl.GetLeafTaxIds();
    // 		if (tax_id_set.empty()) {
    // 			if (dl.IsSetTaxid()) {
    // 				tax_id_set.insert(dl.GetTaxid());
    // 			}
    // 			else {
    // 				tax_id_set.insert(ZERO_TAX_ID);
    // 			}
    // 		}
    //
    // 		string separator = kEmptyStr;
    // 		ITERATE(set<TTaxId>, itr, tax_id_set) {
    // 			if (fields.leaf_node_tax_names == 1) {
    // 				try {
    // 					SSeqDBTaxInfo taxinfo;
    // 					CSeqDB::GetTaxInfo(*itr, taxinfo);
    // 					results[CBlastDeflineUtil::leaf_node_scientific_names] += separator + taxinfo.scientific_name;
    // 					results[CBlastDeflineUtil::leaf_node_common_names] += separator + taxinfo.common_name;
    // 				} catch (const CException&) {
    // 				    results[CBlastDeflineUtil::leaf_node_scientific_names] += separator + NOT_AVAILABLE;
    // 				    results[CBlastDeflineUtil::leaf_node_common_names] += separator + NOT_AVAILABLE;
    // 				}
    // 			}
    // 		    results[CBlastDeflineUtil::leaf_node_tax_ids] += separator + NStr::NumericToString(*itr);
    // 			separator = SEPARATOR;
    // 		}
    // 	}
    //
    // 	if (fields.membership == 1) {
    // 		int membership = 0;
    // 		if(dl.IsSetMemberships()) {
    // 			 ITERATE(CBlast_def_line::TMemberships, memb_int, dl.GetMemberships()) {
    // 				 membership += *memb_int;
    // 			 }
    // 		}
    // 		results[CBlastDeflineUtil::membership] = NStr::NumericToString(membership);
    // 	}
    //
    // 	if (fields.pig == 1) {
    // 		int pig = -1;
    //         if (dl.IsSetOther_info()) {
    //         	ITERATE(CBlast_def_line::TOther_info, itr, dl.GetOther_info()) {
    //         		if (*itr != -1) {
    //         			pig = *itr;
    //         			break;
    // 		        }
    // 		    }
    // 		}
    //         results[CBlastDeflineUtil::pig] = NStr::NumericToString(pig);
    // 	}
    // 	if(fields.links == 1) {
    // 		if (dl.IsSetLinks()) {
    // 			ITERATE(CBlast_def_line::TLinks, links_int, dl.GetLinks()) {
    // 				results[CBlastDeflineUtil::links] += NStr::NumericToString(*links_int) + SEPARATOR;
    // 			}
    // 		}
    // 		else {
    // 			results[CBlastDeflineUtil::links] = NOT_AVAILABLE;
    // 		}
    // 	}
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/blastdb/Blast_def_line.cpp:143-153
    // ```c++
    // CBlast_def_line::TTaxIds CBlast_def_line::GetLeafTaxIds() const
    // {
    //     TTaxIds retval;
    //     if (IsSetLinks()) {
    //         ...
    //         retval.insert(taxids.begin(), taxids.end());
    //     }
    //     return retval;
    // }
    // ```
    /// The fields of `defline` that the specifiers `types` read. The sum of the memberships
    /// is taken in 64 bits (NCBI's `int` would overflow: LOSAT does not reproduce that).
    fn extract(source: &Source, defline: &BlastDefline, types: &[u8]) -> Result<Self, Stop> {
        let wants = |kinds: &[u8]| types.iter().any(|t| kinds.contains(t));
        let mut fields = DeflineFields::default();
        if wants(b"g") {
            fields.gi = match defline.ids.iter().find_map(ParsedSeqId::gi) {
                Some(gi) => gi.to_string(),
                None => NOT_AVAILABLE.to_string(),
            };
        }
        if wants(b"ai") {
            let Some(best) = find_best(&defline.ids, text_score) else {
                // NCBI dereferences the null choice (it crashes on a definition line
                // without Seq-ids).
                return Err(Stop::Refused(format!(
                    "a definition line without Seq-ids is not supported by LOSAT's {}",
                    super::PROGRAM
                )));
            };
            fasta(best, &mut fields.seq_id);
            bare(best, &mut fields.accession);
        }
        if wants(b"t") {
            fields.title = match &defline.title {
                Some(title) => title.clone(),
                None => NOT_AVAILABLE.as_bytes().to_vec(),
            };
        }
        if wants(b"TLBKS") {
            let tax_id = defline.taxid.unwrap_or(0);
            fields.tax_id = tax_id.to_string();
            if wants(b"LBKS") {
                match source.tax_names(tax_id)? {
                    Some(names) => {
                        fields.scientific_name = names.scientific_name;
                        fields.common_name = names.common_name;
                        fields.blast_name = check_name(&names.blast_name);
                        fields.super_kingdom = check_name(&names.super_kingdom);
                    }
                    None => {
                        for field in [
                            &mut fields.scientific_name,
                            &mut fields.common_name,
                            &mut fields.blast_name,
                            &mut fields.super_kingdom,
                        ] {
                            *field = NOT_AVAILABLE.as_bytes().to_vec();
                        }
                    }
                }
            }
        }
        if wants(b"XCN") {
            let mut tax_id_set: BTreeSet<i32> = defline.links.iter().flatten().copied().collect();
            if tax_id_set.is_empty() {
                tax_id_set.insert(defline.taxid.unwrap_or(0));
            }
            let mut separator = "";
            for tax_id in tax_id_set {
                if wants(b"CN") {
                    fields
                        .leaf_node_scientific_names
                        .extend_from_slice(separator.as_bytes());
                    fields
                        .leaf_node_common_names
                        .extend_from_slice(separator.as_bytes());
                    match source.tax_names(tax_id)? {
                        Some(names) => {
                            fields
                                .leaf_node_scientific_names
                                .extend_from_slice(&names.scientific_name);
                            fields
                                .leaf_node_common_names
                                .extend_from_slice(&names.common_name);
                        }
                        None => {
                            fields
                                .leaf_node_scientific_names
                                .extend_from_slice(NOT_AVAILABLE.as_bytes());
                            fields
                                .leaf_node_common_names
                                .extend_from_slice(NOT_AVAILABLE.as_bytes());
                        }
                    }
                }
                fields.leaf_node_tax_ids.push_str(separator);
                fields.leaf_node_tax_ids.push_str(&tax_id.to_string());
                separator = SEPARATOR;
            }
        }
        if wants(b"e") {
            let membership: i64 = defline
                .memberships
                .iter()
                .flatten()
                .map(|&m| i64::from(m))
                .sum();
            fields.membership = membership.to_string();
        }
        if wants(b"P") {
            let pig = defline
                .other_info
                .iter()
                .flatten()
                .copied()
                .find(|&value| value != -1)
                .unwrap_or(-1);
            fields.pig = pig.to_string();
        }
        if wants(b"n") {
            fields.links = match &defline.links {
                Some(links) => links
                    .iter()
                    .map(|link| format!("{link}{SEPARATOR}"))
                    .collect(),
                None => NOT_AVAILABLE.to_string(),
            };
        }
        Ok(fields)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_spec_is_split_as_ncbi_s_loop_does() {
        let spec = FormatSpec::parse("%a%%%l|%%").unwrap();
        assert_eq!(spec.types, b"al");
        assert_eq!(
            spec.separators,
            [b"".to_vec(), b"%".to_vec(), b"|%".to_vec()]
        );
        assert!(spec.get_defline);
        let spec = FormatSpec::parse("%o %l %h").unwrap();
        assert!(!spec.get_defline && spec.need_hash && !spec.need_seq);
        // A `%` at the end reads the terminating NUL.
        assert_eq!(FormatSpec::parse("a%").unwrap().types, [0]);
        assert_eq!(FormatSpec::parse("%").unwrap().separators.len(), 2);
        assert_eq!(
            FormatSpec::parse("no spec").unwrap_err(),
            "Invalid format specification"
        );
        assert_eq!(
            FormatSpec::parse("%%").unwrap_err(),
            "Invalid format specification"
        );
        assert_eq!(
            FormatSpec::parse("").unwrap_err(),
            "Invalid format specification"
        );
    }
}
