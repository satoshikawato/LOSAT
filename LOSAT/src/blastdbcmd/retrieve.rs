//! The retrieval of sequences: NCBI's `x_ProcessSearchRequest` (the configuration of
//! `-outfmt`, `-strand`, `-ctrl_a`, `-mask_sequence_with`, `-range`), its dispatch
//! (`-entry all`, `-entry`, `-entry_batch` with and without `-get_dups`), the ids
//! of a version 5 database (the LMDB table of accessions, the ISAM of GIs), and those LOSAT
//! resolves in a version 4 database without the lookup files: `gnl|BL_ORD_ID|n`, and the
//! answers of a database built without `-parse_seqids`. The lookup of an accession or a GI in
//! a version 4 database built with `-parse_seqids` (ISAM) is line BDC's BDC2c: refused here.

use std::io::{BufWriter, Write};

use crate::seqdb::lookup::IdLookup;
use crate::seqdb::tax_files::NO_ACCESSION_INFO;
use crate::seqdb::taxdb::{LazyTaxDb, TaxNames};
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};
use crate::seqid::accguide::parse_flags::{ANY_LOCAL, DEFAULT, PARTIAL_OK, RAW_TEXT};
use crate::seqid::accguide::SEQID_NOT_SET;
use crate::seqid::label::{content, fasta};
use crate::seqid::model::{Dbtag, ObjectId};
use crate::seqid::parse::{
    parse_ids, seq_id_from_text, string_to_non_negative_int, which_inverse_seq_id,
};
use crate::seqid::rank::{best_rank_score, find_best};
use crate::seqid::{DeflineSet, ParsedSeqId};

use super::args::{self, Args};
use super::fasta::FastaFormatter;
use super::format::FormatSpec;
use super::letters::{protein_fasta_letters, protein_string_letters};
use super::range::{parse_sequence_range_open_end, SeqRange, NO_RANGE};
use super::{Diag, PROGRAM};

/// What ends the retrieval before its last sequence.
#[derive(Debug)]
pub(super) enum Stop {
    /// An exception NCBI's `x_ProcessSearchRequest` catches: its message as a diagnostic
    /// (`Error: [blastdbcmd] <text>`), exit status 1.
    Ncbi(String),
    /// The same, for a message that is not UTF-8 (the unknown specifier read from the NUL
    /// that ends `-outfmt`).
    NcbiBytes(Vec<u8>),
    /// LOSAT's refusal: a diagnostic line, exit status 1.
    Refused(String),
    /// A write to the output failed. NCBI's stream has no exceptions enabled: the run ends
    /// as if nothing had happened, exit status 0.
    Io,
}

impl Stop {
    pub fn io(_: std::io::Error) -> Stop {
        Stop::Io
    }

    fn refused(reason: String) -> Stop {
        Stop::Refused(super::refusal(&reason))
    }
}

impl From<SeqDbError> for Stop {
    /// A reader error inside the request: NCBI's own (`CSeqDBException`) is caught like any
    /// exception there (exit status 1, not the 2 of an error at opening); LOSAT's is a
    /// refusal.
    fn from(error: SeqDbError) -> Stop {
        match error.ncbi_message() {
            Some(text) => Stop::Ncbi(text),
            None => Stop::refused(error.to_string()),
        }
    }
}

/// `CBlastDB_FormatterConfig`, the part LOSAT serves (`m_FiltAlgoId` and `m_FmtAlgoId`
/// never reach a formatter: a database without masks refuses every mask algorithm).
#[derive(Clone, Debug)]
pub(super) struct Config {
    /// `m_SeqRange`.
    pub range: SeqRange,
    /// `m_Strand == eNa_strand_minus` (never for a protein database).
    pub strand_minus: bool,
    /// `m_UseCtrlA`.
    pub ctrl_a: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            range: NO_RANGE,
            strand_minus: false,
            ctrl_a: false,
        }
    }
}

/// The database as the formatters read it.
pub(super) struct Source<'a> {
    pub db: &'a SeqDb,
    /// `m_DbIsProtein`.
    pub protein: bool,
    taxdb: LazyTaxDb,
}

impl<'a> Source<'a> {
    fn new(db: &'a SeqDb, path: &SearchPath) -> Self {
        Source {
            db,
            protein: db.seqtype() == SeqType::Protein,
            taxdb: LazyTaxDb::new(path.clone()),
        }
    }

    /// `GetSeqLength(oid)`.
    pub fn length(&self, oid: u32) -> Result<u32, Stop> {
        if self.protein {
            Ok(self.db.protein_length(oid)?)
        } else {
            let length = self.db.nucleotide_length(oid)?;
            u32::try_from(length)
                .map_err(|_| Stop::refused(format!("sequence {oid} is longer than 32 bits")))
        }
    }

    /// `GetHdr(oid)`: the definition lines with their Seq-ids.
    pub fn deflines(&self, oid: u32) -> Result<DeflineSet, Stop> {
        let set = self.db.defline_set(oid)?;
        if set
            .deflines
            .iter()
            .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
        {
            return Err(Stop::refused(format!(
                "sequence {oid}: a definition line with bytes NCBI reads as '#' (it writes an \
                 error each time it reads them)"
            )));
        }
        Ok(set)
    }

    /// The letters of `GetSequenceAsString(oid[, range])` (`-outfmt %s`, `%h`): iupacaa
    /// with the byte 255 for a gap or stop, or the IUPAC letters of a nucleotide.
    pub fn string_letters(&self, oid: u32, range: Option<(u32, u32)>) -> Result<Vec<u8>, Stop> {
        if self.protein {
            let block = self.db.read_sequence(oid)?;
            Ok(protein_string_letters(slice(block.protein(0), range)))
        } else {
            self.nucleotide_letters(oid, range)
        }
    }

    /// The letters of the FASTA path (`CSeqVector` in `Ncbieaa` or `Iupacna`).
    pub fn fasta_letters(&self, oid: u32, range: Option<(u32, u32)>) -> Result<Vec<u8>, Stop> {
        if self.protein {
            let block = self.db.read_sequence(oid)?;
            Ok(protein_fasta_letters(slice(block.protein(0), range)))
        } else {
            self.nucleotide_letters(oid, range)
        }
    }

    fn nucleotide_letters(&self, oid: u32, range: Option<(u32, u32)>) -> Result<Vec<u8>, Stop> {
        let (from, to) = match range {
            Some((from, to)) => (from as usize, to as usize + 1),
            None => (0, self.db.nucleotide_length(oid)?),
        };
        Ok(self.db.nucleotide_iupac_range(oid, from, to)?)
    }

    /// `CSeqDB::GetTaxInfo(taxid)`: the names of the taxonomy database on the search path,
    /// `None` where NCBI's lookup throws (no database, no such taxid).
    pub fn tax_names(&self, taxid: i32) -> Result<Option<TaxNames>, Stop> {
        let mut messages = Vec::new();
        let names = self.taxdb.names(taxid, &mut messages)?;
        if let Some(message) = messages.first() {
            return Err(Stop::refused(format!(
                "a taxonomy database NCBI opens with a diagnostic ({message})"
            )));
        }
        Ok(names)
    }
}

fn slice(letters: &[u8], range: Option<(u32, u32)>) -> &[u8] {
    match range {
        Some((from, to)) => &letters[from as usize..=to as usize],
        None => letters,
    }
}

/// The formatter of the request (`CBlastDB_Formatter`), or the refusal of one LOSAT does
/// not write, given at the first sequence (NCBI's `%b` and `%d` fail on nothing before).
enum Formatter {
    Fasta(FastaFormatter),
    Fields(FormatSpec),
    Refused(String),
}

impl Formatter {
    fn write(
        &mut self,
        source: &Source,
        oid: u32,
        config: &Config,
        target: Option<&str>,
        out: &mut dyn Write,
    ) -> Result<(), Stop> {
        match self {
            Formatter::Fasta(fasta) => fasta.write(source, oid, config, target, out),
            Formatter::Fields(spec) => spec.write(source, oid, config, target, out),
            Formatter::Refused(reason) => Err(Stop::Refused(reason.clone())),
        }
    }
}

/// `NStr::StringToNumeric<Int8>(text, fConvErr_NoThrow)` with `errno` clear: the number of
/// a text that is one (an optional sign and digits that fit 64 bits).
pub(super) fn numeric_target(text: &str) -> Option<i64> {
    args::string_to_int8_message(text).ok()
}

/// `CSeq_id(text, fParse_PartialOK | fParse_Default)`: the Seq-id of a `-target_only`
/// target that is not a number; the message of NCBI's exception as the request reports it.
pub(super) fn parse_target(text: &str) -> Result<ParsedSeqId, Stop> {
    seq_id_from_text(text.as_bytes(), PARTIAL_OK | DEFAULT, &mut |_| {})
        .map_err(|error| Stop::Ncbi(error.message))
}

/// Runs the retrieval of `args` on the open database `db` (`x_ProcessSearchRequest`), the
/// output to `out`, the diagnostics through `diag`; the exit status.
pub(super) fn process(
    db: &SeqDb,
    args: &Args,
    path: &SearchPath,
    diag: &mut Diag,
    out: &mut dyn Write,
) -> i32 {
    let source = Source::new(db, path);
    let mut sink = BufWriter::new(out);
    let status = match request(args, &source, diag, &mut sink) {
        Ok(false) => 0,
        Ok(true) => 1,
        Err(Stop::Ncbi(message)) => {
            diag.err_post(&message);
            1
        }
        Err(Stop::NcbiBytes(message)) => {
            diag.err_post_bytes(&message);
            1
        }
        Err(Stop::Refused(reason)) => diag.refuse(&reason),
        Err(Stop::Io) => 0,
    };
    let _ = sink.flush();
    status
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:869-896
// ```c++
// int
// CBlastDBCmdApp::x_ProcessSearchRequest()
// {
//    	int err_found = 0;
//     try {
//     	const CArgs& args = GetArgs();
//     	CNcbiOstream& out = args[kArgOutput].AsOutputFile();
//     	string outfmt = x_InitSearchRequest();
//     	/* Special case: full db dump when no range and mask data is specified */
//     	if (m_FASTA) {
//     		CBlastDB_FastaFormatter fasta_fmt(*m_BlastDb, out, args["line_length"].AsInteger(), x_UseLongSeqIds());
//     		err_found = x_ProcessSearchType(fasta_fmt);
//     	}
//     	else if (m_Asn1Bioseq) {
//     		CBlastDB_BioseqFormatter bioseq_fmt(*m_BlastDb, out);
//     		err_found = x_ProcessSearchType(bioseq_fmt);
//     	}
//     	else {
//     		CBlastDB_SeqFormatter seq_fmt(outfmt, *m_BlastDb, out);
//     		err_found = x_ProcessSearchType(seq_fmt);
//     	}
//     }
//     catch (const CException& e) {
//     	ERR_POST(Error << e.GetMsg());
//         err_found = 1;
//     } catch (...) {
//         ERR_POST(Error << "Failed to retrieve requested item");
//         err_found = 1;
//     }
// 	return err_found;
// }
// ```
/// The request: whether an entry was not found (`err_found`), or what stopped it.
fn request(
    args: &Args,
    source: &Source,
    diag: &mut Diag,
    out: &mut dyn Write,
) -> Result<bool, Stop> {
    let (mut config, kind) = init_search_request(args, source, diag)?;
    let line_length = args
        .value("line_length")
        .and_then(|value| value.parse::<u32>().ok())
        .unwrap_or(80);
    let mut formatter = match kind {
        OutputKind::Fasta => {
            Formatter::Fasta(FastaFormatter::new(line_length, args.given("long_seqids")))
        }
        OutputKind::Bioseq => Formatter::Refused(super::refusal(
            "-outfmt %b (the Bioseq as ASN.1 text) is not served",
        )),
        OutputKind::Fields(spec) => match FormatSpec::parse(&spec) {
            Ok(spec) if spec.asks_asn_text() => Formatter::Refused(super::refusal(
                "-outfmt %d (the definition lines as ASN.1 text) is not served",
            )),
            Ok(spec) => Formatter::Fields(spec),
            Err(message) => return Err(Stop::Ncbi(message)),
        },
    };
    let ids = DbIds::of(source, args)?;
    let target_only = args.given("target_only");
    let mut run = Run {
        source,
        ids: &ids,
        formatter: &mut formatter,
        target_only,
        diag,
        out,
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:819-846
    // ```c++
    // 	if (args["entry"].HasValue() && args["entry"].AsString() == "all") {
    // 		fmt.DumpAll(m_Config);
    // 	}
    // 	else if (args["entry_batch"].HasValue()) {
    // 		if(m_GetDuplicates) {
    // 			return x_ProcessBatchEntry(fmt);
    // 		}
    // 		else {
    // 			return x_ProcessBatchEntry_NoDup(fmt);
    // 		}
    // 	}
    // 	else if (args["entry"].HasValue() || args["ipg"].HasValue()) {
    // 		return x_ProcessEntry(fmt);
    // 	}
    // 	...
    // 	else {
    // 		NCBI_THROW(CInputException, eInvalidInput,
    // 		       	   "Must specify query type: one of 'entry', 'entry_batch', or 'pig'");
    // 	}
    // ```
    match (args.value("entry"), args.value("entry_batch")) {
        (Some("all"), _) => {
            run.dump_all(&config)?;
            Ok(false)
        }
        (_, Some(file)) => {
            if args.given("get_dups") {
                run.batch_entry(file, &mut config)
            } else {
                run.batch_entry_no_dup(file, &mut config)
            }
        }
        (Some(entry), None) => run.entry(entry, &config),
        (None, None) => Err(Stop::Ncbi(
            "Must specify query type: one of 'entry', 'entry_batch', or 'pig'".to_string(),
        )),
    }
}

/// What `x_InitSearchRequest` chooses: `m_FASTA`, `m_Asn1Bioseq`, or the format text of the
/// `%` formatter.
enum OutputKind {
    Fasta,
    Bioseq,
    Fields(String),
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:736-816
// ```c++
// string
// CBlastDBCmdApp::x_InitSearchRequest()
// {
//    	const CArgs& args = GetArgs();
//     m_GetDuplicates = args["get_dups"];
//     m_TargetOnly = args["target_only"];
//
//     string outfmt = kEmptyStr;
//     if (args["outfmt"].HasValue()) {
//     	outfmt = args["outfmt"].AsString();
//         m_FASTA = false;
//         m_Asn1Bioseq = false;
//
//         if ((outfmt.find("%f") != string::npos &&
//            	(outfmt.find("%b") != string::npos || outfmt.find("%d") != string::npos)) ||
//             (outfmt.find("%b") != string::npos && outfmt.find("%d") != string::npos)) {
//            	NCBI_THROW(CInputException, eInvalidInput,
//                     	"The %f, %b, %d output format options cannot be specified together.");
//         }
//
//         if (outfmt.find("%b") != string::npos) {
//            	outfmt = "%b";
//            	m_Asn1Bioseq = true;
//         }
//
//         // If "%f" is found within outfmt, discard everything else
//         if (outfmt.find("%f") != string::npos) {
//            	outfmt = "%f";
//            	m_FASTA = true;
//         }
//
//         if (outfmt.find("%d") != string::npos) {
//            	outfmt = "%d";
//         }
//
//         if (outfmt.find("%m") != string::npos) {
//            	int algo_id = 0;
//            	size_t i = outfmt.find("%m") + 2;
//            	bool found = false;
//            	while (i < outfmt.size() && outfmt[i] >= '0' && outfmt[i] <= '9') {
//            		algo_id = algo_id * 10 + (outfmt[i] - '0');
//                	outfmt.erase(i, 1);
//                	found = true;
//             }
//             if (!found) {
//                	NCBI_THROW(CInputException, eInvalidInput,
//                        	   "The option '-outfmt %m' is not followed by a masking algo ID.");
//             }
//             m_Config.m_FmtAlgoId = algo_id;
//     		if(!s_IsMaskAlgoIdValid(*m_BlastDb, m_Config.m_FmtAlgoId)) {
//     			NCBI_THROW(CInvalidDataException, eInvalidInput,
//     				                   "Invalid filtering algorithm ID for outfmt %m.");
//     		}
//         }
//     }
//
//     if (args["strand"].HasValue() && !m_DbIsProtein) {
//     	if (args["strand"].AsString() == "plus") {
//                 m_Config.m_Strand = eNa_strand_plus;
//             } else if (args["strand"].AsString() == "minus") {
//                 m_Config.m_Strand = eNa_strand_minus;
//             } else {
//             	NCBI_THROW(CInputException, eInvalidInput,
//             	           "Both strands is not supported");
//             }
//     }
//     m_Config.m_UseCtrlA = args["ctrl_a"];
//     if (args["mask_sequence_with"].HasValue()) {
//     	m_Config.m_FiltAlgoId = -1;
//         m_Config.m_FiltAlgoId = NStr::StringToInt(args["mask_sequence_with"].AsString(), NStr::fConvErr_NoThrow);
//         if(errno) {
//         	m_Config.m_FiltAlgoId = m_BlastDb->GetMaskAlgorithmId(args["mask_sequence_with"].AsString());
//         }
//    		if(!s_IsMaskAlgoIdValid(*m_BlastDb, m_Config.m_FiltAlgoId)){
//    			NCBI_THROW(CInvalidDataException, eInvalidInput,
//    		               "Invalid filtering algorithm ID for mask_sequence_with.");
//    		}
//     }
//     if (args["range"].HasValue()) {
//     	m_Config.m_SeqRange = ParseSequenceRangeOpenEnd(args["range"].AsString());
//     }
//      return outfmt;
// }
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:259-270
// ```c++
// bool s_IsMaskAlgoIdValid(CSeqDB & blastdb, int id)
// {
// 	if (id >= 0) {
// 	    vector<int> algo_id(1, id);
// 	    vector<int> invalid_algo_ids = blastdb.ValidateMaskAlgorithms(algo_id);
// 	    if ( !invalid_algo_ids.empty()) {
// 	    	ERR_POST(Error << "Invalid filtering algorithm ID: " << NStr::IntToString(id));
// 	    	return false;
// 	    }
// 	}
// 	return true;
// }
// ```
/// `x_InitSearchRequest`: the configuration and the kind of output. Every mask algorithm id
/// is invalid: LOSAT's reader refuses the databases that have masks (the column files), so
/// the database here has none, and NCBI's answers for it are reproduced.
fn init_search_request(
    args: &Args,
    source: &Source,
    diag: &mut Diag,
) -> Result<(Config, OutputKind), Stop> {
    let mut config = Config::default();
    let mut kind = OutputKind::Fasta;
    if let Some(given) = args.value("outfmt") {
        let mut outfmt = given.to_string();
        let has = |text: &str, spec: &str| text.contains(spec);
        if (has(&outfmt, "%f") && (has(&outfmt, "%b") || has(&outfmt, "%d")))
            || (has(&outfmt, "%b") && has(&outfmt, "%d"))
        {
            return Err(Stop::Ncbi(
                "The %f, %b, %d output format options cannot be specified together.".to_string(),
            ));
        }
        kind = OutputKind::Fields(String::new());
        if has(&outfmt, "%b") {
            outfmt = "%b".to_string();
            kind = OutputKind::Bioseq;
        }
        if has(&outfmt, "%f") {
            outfmt = "%f".to_string();
            kind = OutputKind::Fasta;
        }
        if has(&outfmt, "%d") {
            outfmt = "%d".to_string();
        }
        if let Some(at) = outfmt.find("%m") {
            let digits_start = at + 2;
            let digits_end = outfmt[digits_start..]
                .find(|c: char| !c.is_ascii_digit())
                .map_or(outfmt.len(), |len| digits_start + len);
            if digits_end == digits_start {
                return Err(Stop::Ncbi(
                    "The option '-outfmt %m' is not followed by a masking algo ID.".to_string(),
                ));
            }
            let algo_id = outfmt[digits_start..digits_end]
                .bytes()
                .fold(0i64, |id, digit| {
                    id.wrapping_mul(10).wrapping_add(i64::from(digit - b'0'))
                });
            outfmt.replace_range(digits_start..digits_end, "");
            if algo_id >= 0 {
                diag.err_post(&format!("Invalid filtering algorithm ID: {algo_id}"));
                return Err(Stop::Ncbi(
                    "Invalid filtering algorithm ID for outfmt %m.".to_string(),
                ));
            }
        }
        if let OutputKind::Fields(text) = &mut kind {
            *text = outfmt;
        }
    }
    if !source.protein && args.value("strand") == Some("minus") {
        config.strand_minus = true;
    }
    config.ctrl_a = args.given("ctrl_a");
    if let Some(value) = args.value("mask_sequence_with") {
        match args::string_to_int8_message(value)
            .ok()
            .and_then(|id| i32::try_from(id).ok())
        {
            Some(id) if id >= 0 => {
                diag.err_post(&format!("Invalid filtering algorithm ID: {id}"));
                return Err(Stop::Ncbi(
                    "Invalid filtering algorithm ID for mask_sequence_with.".to_string(),
                ));
            }
            Some(id) => {
                return Err(Stop::refused(format!(
                    "-mask_sequence_with {id} (a negative algorithm id, which NCBI accepts \
                     and then looks up) is not served"
                )));
            }
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2219-2224
            // ```c++
            //     if (algo_id_map.find(algo_name) == algo_id_map.end()) {
            //         NCBI_THROW(CSeqDBException, eArgErr, "Cannot find string algorithm id in algorithm map.");
            //     }
            // ```
            None => {
                return Err(Stop::Ncbi(
                    "Cannot find string algorithm id in algorithm map.".to_string(),
                ));
            }
        }
    }
    if let Some(value) = args.value("range") {
        if value.is_empty() {
            // `tokens.front()` of an empty vector: NCBI crashed (probe D61).
            return Err(Stop::refused(
                "-range \"\" (NCBI reads past an empty list of range elements) is not served"
                    .to_string(),
            ));
        }
        config.range = parse_sequence_range_open_end(value).map_err(Stop::Ncbi)?;
    }
    Ok((config, kind))
}

/// What the database's files say about looking ids up: which lookups LOSAT answers.
struct DbIds<'a> {
    /// The tables of a version 5 database, when the request looks ids up (the volume lists
    /// LOSAT refuses do not stop `-entry all`).
    lookup: Option<IdLookup<'a>>,
    /// Version 4: a volume has the numeric ISAM of GIs (`.pni`/`.nni`).
    gi_isam: bool,
    /// Version 4: a volume has the string ISAM (`.psi`/`.nsi`).
    str_isam: bool,
    num_seqs: u32,
}

impl<'a> DbIds<'a> {
    fn of(source: &Source<'a>, args: &Args) -> Result<Self, Stop> {
        let db = source.db;
        let mut ids = DbIds {
            lookup: None,
            gi_isam: false,
            str_isam: false,
            num_seqs: db.num_seqs(),
        };
        if db.is_version5() {
            let looks_up = args.given("entry_batch")
                || args.value("entry").is_some_and(|entry| entry != "all");
            if looks_up {
                ids.lookup = Some(db.id_lookup()?);
            }
        } else {
            let letter = if source.protein { "p" } else { "n" };
            for volume in db.volumes() {
                let file = |extension: &str| {
                    let mut name = volume.base().as_os_str().to_os_string();
                    name.push(format!(".{letter}{extension}"));
                    std::fs::metadata(name).is_ok_and(|meta| meta.is_file())
                };
                ids.gi_isam |= file("ni");
                ids.str_isam |= file("si");
            }
        }
        Ok(ids)
    }
}

/// The answer of `x_GetOids` for one id.
enum Lookup {
    Oids(Vec<u32>),
    /// `Entry not found: <acc>` (the id as NCBI's lookup saw it).
    NotFound(String),
    /// An exception of the lookup with NCBI's text, or LOSAT's refusal.
    Stop(Stop),
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:152-180
// ```c++
// string s_PreProcessAccessionsForDBv5(const string & id)
// {
// 	string rv = id;
// 	if ((id.find('|') != NPOS) || (id.find('_') != NPOS)) {
//
// 		CRef<CSeq_id> seqid;
// 		try {
// 			seqid = new CSeq_id(id, CSeq_id::fParse_RawText | CSeq_id::fParse_AnyLocal | CSeq_id::fParse_PartialOK);
// 		}
// 		catch(...) {
// 		}
//
// 		if(seqid.NotEmpty()) {
// 			if(seqid->IsPir() || seqid->IsPrf()) {
// 				return seqid->AsFastaString();
// 			}
// 			else if (seqid->IsPdb()) {
// 				string tmp = seqid->GetSeqIdString();
// 				rv = tmp.substr(0,4);
// 				rv += tmp.substr(4);
// 				return (rv);
// 			}
// 			return seqid->GetSeqIdString(true);
// 		}
// 	}
//
// 	return NStr::ToUpper(rv);
//
// }
// ```
/// The key a version 5 database is asked for: the label of an id with `|` or `_`, else the
/// text in upper case.
fn preprocess_accession_for_v5(id: &str) -> String {
    if id.contains('|') || id.contains('_') {
        if let Ok(seqid) = seq_id_from_text(
            id.as_bytes(),
            RAW_TEXT | ANY_LOCAL | PARTIAL_OK,
            &mut |_| {},
        ) {
            let mut label = Vec::new();
            match &seqid {
                ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_) => fasta(&seqid, &mut label),
                ParsedSeqId::Pdb(_) => content(&seqid, false, &mut label),
                _ => content(&seqid, true, &mut label),
            }
            return String::from_utf8_lossy(&label).into_owned();
        }
    }
    id.to_ascii_uppercase()
}

/// What `SeqDB_SimplifyAccession` makes of an id of a version 4 database, as far as the
/// lookup that follows it: an ordinal id (no file read), the exception of a `BL_ORD_ID`
/// whose tag is a string, or a key of the ISAM files.
enum Simplified {
    Oid(i32),
    StringTag,
    IsamKey,
    ParseError(String),
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2535-2548
// ```c++
// ESeqDBIdType SeqDB_SimplifyAccession(const string & acc,
//                                      Int8         & num_id,
//                                      string       & str_id,
//                                      bool         & simpler)
// {
//     ESeqDBIdType result = eStringId;
//     num_id = (Uint4)-1;
//
//     vector< CRef< CSeq_id > > seqid_set;
//
//     if (s_SeqDB_ParseSeqIDs(acc, seqid_set)) {
//         // Something like SeqIdFindBest()
//         CRef<CSeq_id> bestid =
//             FindBestChoice(seqid_set, CSeq_id::BestRank);
//
//         result = SeqDB_SimplifySeqid(*bestid, & acc, num_id, str_id, simpler);
//     } else {
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2286-2296
// ```c++
//     case CSeq_id::e_General:
//         {
//             const CDbtag & dbt = bestid.GetGeneral();
//
//             if (dbt.CanGetDb()) {
//                 if (dbt.GetDb() == "BL_ORD_ID") {
//                     simpler = true;
//                     num_id = dbt.GetTag().GetId();
//                     result = eOID;
//                     break;
//                 }
// ```
/// `s_SeqDB_ParseSeqIDs` reads ids only from a text that starts with a Seq-id tag before
/// its first `|`; the best by `BestRank` goes to `SeqDB_SimplifySeqid`.
fn simplify_v4(acc: &str) -> Simplified {
    let Some(bar) = acc.find('|') else {
        return Simplified::IsamKey;
    };
    if which_inverse_seq_id(acc[..bar].as_bytes()) == SEQID_NOT_SET {
        return Simplified::IsamKey;
    }
    let ids: Vec<ParsedSeqId> = match parse_ids(acc.as_bytes(), DEFAULT, &mut |_| {}) {
        Ok(ids) => ids.into_iter().flatten().collect(),
        Err(error) => return Simplified::ParseError(error.message),
    };
    match find_best(&ids, best_rank_score) {
        Some(ParsedSeqId::General(Dbtag { db, tag })) if db == b"BL_ORD_ID" => match tag {
            ObjectId::Id(n) => Simplified::Oid(*n),
            // `CObject_id::GetId()` of a string tag: `CInvalidChoiceSelection`, whose message
            // the request reports (probe D11, D51-D55: the tags `0`, `-1`, `abc`, `03`).
            ObjectId::Str(_) => Simplified::StringTag,
        },
        _ => Simplified::IsamKey,
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:183-210
// ```c++
// bool
// CBlastDBCmdApp::x_GetOids(const string & id, vector<int> & oids)
// {
// 	string acc = id;
// 	if(m_BlastDb->GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5) {
// 		acc = s_PreProcessAccessionsForDBv5(id);
// 	}
// 	TGi num_id = NStr::StringToNumeric<TGi>(acc, NStr::fConvErr_NoThrow);
// 	if(!errno) {
// 		int gi_oid = -1;
// 		m_BlastDb->GiToOidwFilterCheck(num_id, gi_oid);
// 		if(gi_oid < 0) {
// 			m_BlastDb->AccessionToOids(acc, oids);
// 		}
// 		else {
// 			oids.push_back(gi_oid);
// 		}
//
// 	}
// 	else {
// 		m_BlastDb->AccessionToOids(acc, oids);
// 	}
// 	if(oids.empty()) {
// 		ERR_POST(Error <<  "Entry not found: " << acc);
// 		return false;
// 	}
// 	return true;
// }
// ```
/// `x_GetOids(id)` as LOSAT answers it. A version 5 database looks the key up in its tables
/// (`get_oids_v5`); `gnl|BL_ORD_ID|n` of a version 4 database is the ordinal id `n` with no
/// file read (line BDC DECISIONS 12: the one sequence `n` of the whole database, where NCBI
/// adds `n` to the start of every volume); the other ids of a version 4 database are not
/// found when it has no lookup files. The lookups in its files are refused.
fn get_oids(ids: &DbIds, id: &str) -> Lookup {
    if let Some(lookup) = &ids.lookup {
        return get_oids_v5(lookup, id);
    }
    let refused = |what: &str| {
        Lookup::Stop(Stop::Refused(super::refusal(&format!(
            "the lookup of {what} in a database built with -parse_seqids (-entry {id}) is not \
             served"
        ))))
    };
    if numeric_target(id).is_some() && ids.gi_isam {
        return refused("a GI");
    }
    match simplify_v4(id) {
        Simplified::Oid(n) => {
            if u32::try_from(n).is_ok_and(|n| n < ids.num_seqs) {
                Lookup::Oids(vec![n as u32])
            } else {
                Lookup::NotFound(id.to_string())
            }
        }
        Simplified::StringTag => Lookup::Stop(Stop::Ncbi(
            "CObject_id::GetId(): Invalid choice selection: NCBI-General::Object-id.str"
                .to_string(),
        )),
        Simplified::ParseError(message) => Lookup::Stop(Stop::Ncbi(message)),
        Simplified::IsamKey => {
            if ids.str_isam {
                refused("an accession")
            } else {
                Lookup::NotFound(id.to_string())
            }
        }
    }
}

/// `x_GetOids(id)` of a version 5 database: a GI first when the key is a number, then the
/// accession table.
fn get_oids_v5(lookup: &IdLookup, id: &str) -> Lookup {
    let acc = preprocess_accession_for_v5(id);
    if let Some(gi) = numeric_target(&acc) {
        match lookup.gi_to_oid(gi) {
            Ok(Some(oid)) => return Lookup::Oids(vec![oid]),
            Ok(None) => {}
            Err(error) => return Lookup::Stop(error.into()),
        }
    }
    match lookup.accession_to_oids(acc.as_bytes()) {
        Ok(oids) if oids.is_empty() => Lookup::NotFound(acc),
        Ok(oids) => Lookup::Oids(oids),
        Err(error) => Lookup::Stop(error.into()),
    }
}

/// The lines of an `-entry_batch` file: the id and the rest of the line
/// (`NStr::SplitInTwo(line, " \t", id, format, fSplit_MergeDelimiters | fSplit_Truncate)`),
/// empty lines and lines without an id left out. The line's `\r` stays, as in NCBI's reading.
fn read_batch(file: &str) -> Result<Vec<(String, String)>, Stop> {
    let bytes = std::fs::read(file).map_err(|_| {
        Stop::Ncbi(args::arg_message(
            "entry_batch",
            "File is not accessible",
            file,
        ))
    })?;
    let text = String::from_utf8(bytes)
        .map_err(|_| Stop::refused(format!("{file}: a batch file that is not UTF-8 text")))?;
    let delimiter = |c: char| c == ' ' || c == '\t';
    let mut lines = Vec::new();
    for line in text.split('\n') {
        let trimmed = line.trim_matches(delimiter);
        if trimmed.is_empty() {
            continue;
        }
        let (id, format) = match trimmed.find(delimiter) {
            Some(at) => (&trimmed[..at], trimmed[at..].trim_start_matches(delimiter)),
            None => (trimmed, ""),
        };
        lines.push((id.to_string(), format.to_string()));
    }
    Ok(lines)
}

/// The retrieval in progress.
struct Run<'a, 'b> {
    source: &'a Source<'a>,
    ids: &'a DbIds<'a>,
    formatter: &'a mut Formatter,
    target_only: bool,
    diag: &'a mut Diag<'b>,
    out: &'a mut dyn Write,
}

impl Run<'_, '_> {
    fn write(&mut self, oid: u32, config: &Config, target: Option<&str>) -> Result<(), Stop> {
        self.formatter
            .write(self.source, oid, config, target, self.out)
    }

    /// Refuses before any output when an id of `ids` would need a lookup LOSAT does not do
    /// (in a version 4 database).
    fn refuse_lookups(&self, ids: impl Iterator<Item = impl AsRef<str>>) -> Result<(), Stop> {
        if self.ids.lookup.is_some() {
            return Ok(());
        }
        for id in ids {
            if let Lookup::Stop(Stop::Refused(reason)) = get_oids(self.ids, id.as_ref()) {
                return Err(Stop::Refused(reason));
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:411-416
    // ```c++
    // void CBlastDB_SeqFormatter::DumpAll(const CBlastDB_FormatterConfig & config)
    // {
    // 	for (int i = 0; m_BlastDb.CheckOrFindOID(i); i++) {
    // 		Write(i, config);
    // 	}
    // }
    // ```
    /// `-entry all`: every sequence in the order of its ordinal id.
    fn dump_all(&mut self, config: &Config) -> Result<(), Stop> {
        for oid in 0..self.ids.num_seqs {
            self.write(oid, config, None)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:213-257
    // ```c++
    //     } else if (args["entry"].HasValue()) {
    //     	static const string kDelim(",");
    //     	const string& entry = args["entry"].AsString();
    //
    //     	vector<string> queries;
    //        	if (entry.find(kDelim[0]) != string::npos) {
    //            	NStr::Split(entry, kDelim, queries);
    //        	} else {
    //        		queries.resize(1);
    //        		queries[0]  = entry;
    //      	}
    //    		for(unsigned int i=0; i < queries.size(); i++) {
    //      		vector<CSeqDB::TOID> oids;
    //      		if(x_GetOids(queries[i], oids)) {
    //    				for(unsigned int j=0; j < oids.size(); j++) {
    //    					if(m_TargetOnly) {
    //    						fmt.Write(oids[j], m_Config, queries[i]);
    //    					}
    //    					else {
    //    						fmt.Write(oids[j], m_Config);
    //    					}
    //      			}
    //      		}
    //      		else {
    //      			err_found ++;
    //      		}
    //      	}
    //    		if(err_found == queries.size()) {
    //    			NCBI_THROW(CInputException, eInvalidInput,
    //    		               "Entry or entries not found in BLAST database");
    //    		}
    //     }
    //    	return (err_found) ? 1:0;
    // ```
    /// `-entry a,b,...`: whether an entry was not found.
    fn entry(&mut self, entry: &str, config: &Config) -> Result<bool, Stop> {
        let queries: Vec<&str> = if entry.contains(',') {
            entry.split(',').collect()
        } else {
            vec![entry]
        };
        self.refuse_lookups(queries.iter())?;
        let mut err_found = 0;
        for query in &queries {
            match get_oids(self.ids, query) {
                Lookup::Oids(oids) => {
                    for oid in oids {
                        self.write(oid, config, self.target_only.then_some(query))?;
                    }
                }
                Lookup::NotFound(acc) => {
                    self.diag.err_post(&format!("Entry not found: {acc}"));
                    err_found += 1;
                }
                Lookup::Stop(stop) => return Err(stop),
            }
        }
        if err_found == queries.len() {
            return Err(Stop::Ncbi(
                "Entry or entries not found in BLAST database".to_string(),
            ));
        }
        Ok(err_found > 0)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:375-438
    // ```c++
    // int
    // CBlastDBCmdApp::x_ProcessBatchEntry_NoDup(CBlastDB_Formatter & fmt)
    // {
    // 	int err_found = 0;
    //    	const CArgs& args = GetArgs();
    //     CNcbiIstream& input = args["entry_batch"].AsInputFile();
    //     vector<string> ids, formats;
    //     vector<CSeqDB::TOID> oids;
    //     while (input) {
    //         ...
    //     }
    //
    //     if(m_BlastDb->GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5) {
    //     	for(unsigned int i=0; i < ids.size(); i++) {
    //     		ids[i] = s_PreProcessAccessionsForDBv5(ids[i]);
    //     	}
    //     }
    //     try {
    //     m_BlastDb->AccessionsToOids(ids, oids);
    //     }
    //     catch (CSeqDBException & e) {
    //     	if (e.GetMsg().find("DB contains no accession info") == NPOS){
    //     		NCBI_RETHROW_SAME(e, e.GetMsg());
    //     	}
    //     }
    //     for(unsigned i=0; i < ids.size(); i++) {
    //     	if(oids[i] == kSeqDBEntryNotFound) {
    //     		TGi num_id = NStr::StringToNumeric<TGi>(ids[i], NStr::fConvErr_NoThrow);
    //     		if(!errno) {
    //     			int gi_oid = -1;
    //     			m_BlastDb->GiToOidwFilterCheck(num_id, gi_oid);
    //     			if(gi_oid >= 0) {
    //     				oids[i] = gi_oid;
    //     			}
    //     		}
    //     		if(oids[i] == kSeqDBEntryNotFound) {
    //     			err_found ++;
    //     			ERR_POST (Error << "Skipped " << ids[i]);
    //     			continue;
    //     		}
    //     	}
    //     	if(x_ModifyConfigForBatchEntry(formats[i]))  {
    //     		err_found ++;
    //     		ERR_POST (Error << "Skipped " << ids[i]);
    //     		continue;
    //     	}
    //     	if(m_TargetOnly) {
    //     		fmt.Write(oids[i], m_Config, ids[i]);
    //     	}
    //     	else {
    //     	   	fmt.Write(oids[i], m_Config);
    //     	}
    //     }
    //     return (err_found) ? 1 : 0;
    // }
    // ```
    /// `-entry_batch` without `-get_dups`: one sequence per id, the ids resolved all at once
    /// first (the exception of a version 5 database without `acc2oid` is swallowed there:
    /// the ids keep what the LMDB files before it answered), a GI after that.
    fn batch_entry_no_dup(&mut self, file: &str, config: &mut Config) -> Result<bool, Stop> {
        let lines = read_batch(file)?;
        let ids: Vec<String> = lines
            .iter()
            .map(|(id, _)| {
                if self.ids.lookup.is_some() {
                    preprocess_accession_for_v5(id)
                } else {
                    id.clone()
                }
            })
            .collect();
        self.refuse_lookups(ids.iter())?;
        let mut oids = Vec::with_capacity(ids.len());
        if let Some(lookup) = &self.ids.lookup {
            let keys: Vec<Vec<u8>> = ids.iter().map(|id| id.as_bytes().to_vec()).collect();
            if let Err(error) = lookup.accessions_to_oids(&keys, &mut oids) {
                let swallowed = error
                    .ncbi_message()
                    .is_some_and(|message| message.contains(NO_ACCESSION_INFO));
                if !swallowed {
                    return Err(error.into());
                }
            }
        } else {
            for id in &ids {
                oids.push(match get_oids(self.ids, id) {
                    Lookup::Oids(found) => Some(found[0]),
                    Lookup::NotFound(_) => None,
                    Lookup::Stop(stop) => return Err(stop),
                });
            }
        }
        let mut err_found = false;
        for (i, id) in ids.iter().enumerate() {
            let mut oid = oids[i];
            if let (None, Some(lookup), Some(gi)) = (oid, &self.ids.lookup, numeric_target(id)) {
                oid = lookup.gi_to_oid(gi)?;
            }
            let Some(oid) = oid else {
                err_found = true;
                self.diag.err_post(&format!("Skipped {id}"));
                continue;
            };
            if self.modify_config_for_batch_entry(&lines[i].1, config) != 0 {
                err_found = true;
                self.diag.err_post(&format!("Skipped {id}"));
                continue;
            }
            self.write(oid, config, self.target_only.then_some(id))?;
        }
        Ok(err_found)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:441-485
    // ```c++
    // int
    // CBlastDBCmdApp::x_ProcessBatchEntry(CBlastDB_Formatter & fmt)
    // {
    // 	int err_found = 0;
    //    	const CArgs& args = GetArgs();
    //     CNcbiIstream& input = args["entry_batch"].AsInputFile();
    //
    //     while (input) {
    //         string line;
    //         NcbiGetlineEOL(input, line);
    //         if ( !line.empty() ) {
    //         	string id, format;
    //         	NStr::SplitInTwo(line, " \t", id, format, NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
    //         	if(id.empty()) {
    //         		continue;
    //         	}
    //         	if(x_ModifyConfigForBatchEntry(format))  {
    //         		err_found ++;
    //         		ERR_POST (Error << "Skipped " << id);
    //         		continue;
    //         	}
    //    			vector<int> oids;
    //          	if(!x_GetOids(id, oids)) {
    //          		err_found ++;
    //         		ERR_POST (Error << "Skipped " << id);
    //         		continue;
    //          	}
    //
    //          	if (m_GetDuplicates) {
    //          		for(unsigned int j=0; j < oids.size(); j++) {
    //          			fmt.Write(oids[j], m_Config);
    //          		}
    //            	}
    //          	...
    //         }
    //     }
    //     return (err_found) ? 1 : 0;
    // }
    // ```
    /// `-entry_batch -get_dups`: line by line, every sequence of each id.
    fn batch_entry(&mut self, file: &str, config: &mut Config) -> Result<bool, Stop> {
        let lines = read_batch(file)?;
        self.refuse_lookups(lines.iter().map(|(id, _)| id))?;
        let mut err_found = false;
        for (id, format) in &lines {
            if self.modify_config_for_batch_entry(format, config) != 0 {
                err_found = true;
                self.diag.err_post(&format!("Skipped {id}"));
                continue;
            }
            match get_oids(self.ids, id) {
                Lookup::Oids(oids) => {
                    for oid in oids {
                        self.write(oid, config, None)?;
                    }
                }
                Lookup::NotFound(acc) => {
                    self.diag.err_post(&format!("Entry not found: {acc}"));
                    err_found = true;
                    self.diag.err_post(&format!("Skipped {id}"));
                }
                Lookup::Stop(stop) => return Err(stop),
            }
        }
        Ok(err_found)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:271-300
    // ```c++
    // int CBlastDBCmdApp::x_ModifyConfigForBatchEntry(const string & format)
    // {
    // 	int status = 0;
    // 	if (!m_DbIsProtein) {
    // 		m_Config.m_Strand = eNa_strand_plus;
    // 	}
    //    	m_Config.m_SeqRange = TSeqRange::GetEmpty();
    //    	m_Config.m_FiltAlgoId = -1;
    //    	if(!format.empty()) {
    //    		vector<string> tmp;
    //    		NStr::Split(format, " \t", tmp, NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
    //    		for(unsigned int i=0; i < tmp.size(); i++) {
    //    			if(tmp[i].find('-')!= string::npos) {
    //    				try {
    //    					m_Config.m_SeqRange = ParseSequenceRangeOpenEnd(tmp[i]);
    //    				} catch (...) {
    //    				}
    //    			}
    //    			else if (!m_DbIsProtein && NStr::EqualNocase(tmp[i].c_str(), "minus")) {
    //    				m_Config.m_Strand = eNa_strand_minus;
    //    			}
    //    			else {
    //    				m_Config.m_FiltAlgoId = NStr::StringToNonNegativeInt(tmp[i]);
    //    				if(!s_IsMaskAlgoIdValid(*m_BlastDb, m_Config.m_FiltAlgoId)){
    //    					status = 1;
    //    				}
    //    			}
    //    		}
    //    	}
    //    	return status;
    // }
    // ```
    /// The range, strand and mask algorithm of a batch line (`id [range] [minus] [algo]`):
    /// 1 when the algorithm id is invalid (every one is, in a database without masks; a word
    /// that is not a number is the id -1, which NCBI does not check).
    fn modify_config_for_batch_entry(&mut self, format: &str, config: &mut Config) -> i32 {
        let mut status = 0;
        if !self.source.protein {
            config.strand_minus = false;
        }
        config.range = NO_RANGE;
        for token in format.split([' ', '\t']).filter(|token| !token.is_empty()) {
            if token.contains('-') {
                if let Ok(range) = parse_sequence_range_open_end(token) {
                    config.range = range;
                }
            } else if !self.source.protein && token.eq_ignore_ascii_case("minus") {
                config.strand_minus = true;
            } else {
                let algo_id = string_to_non_negative_int(token.as_bytes());
                if algo_id >= 0 {
                    self.diag
                        .err_post(&format!("Invalid filtering algorithm ID: {algo_id}"));
                    status = 1;
                }
            }
        }
        status
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version5_keys_are_preprocessed_as_ncbi_s() {
        let key = preprocess_accession_for_v5;
        assert_eq!(key("p69905"), "P69905");
        assert_eq!(key("np_000002.1"), "NP_000002.1");
        assert_eq!(key("sp|P99999.2|TEST_HUMAN"), "P99999.2");
        assert_eq!(key("pdb|1ABC|A"), "1ABC_A");
        assert_eq!(key("28zx_m"), "28zx_m");
        assert_eq!(key("lcl|myprot9"), "myprot9");
        assert_eq!(key("gnl|mydb|prot10"), "mydb:prot10");
        assert_eq!(key("gnl|BL_ORD_ID|3"), "BL_ORD_ID:3");
        assert_eq!(key("pir||s00007"), "pir||s00007");
        assert_eq!(key("gi|123456"), "123456");
        assert_eq!(key("mydb:prot10"), "MYDB:PROT10");
    }

    #[test]
    fn ordinal_ids_of_version4_are_read_without_files() {
        let ids = DbIds {
            lookup: None,
            gi_isam: false,
            str_isam: false,
            num_seqs: 20,
        };
        let oids = |id: &str| match get_oids(&ids, id) {
            Lookup::Oids(oids) => format!("oids {oids:?}"),
            Lookup::NotFound(acc) => format!("not found {acc}"),
            Lookup::Stop(Stop::Ncbi(message)) => format!("error {message}"),
            Lookup::Stop(stop) => format!("{stop:?}"),
        };
        assert_eq!(oids("gnl|BL_ORD_ID|3"), "oids [3]");
        assert_eq!(oids("gnl|BL_ORD_ID|3|"), "oids [3]");
        assert_eq!(oids("gnl|BL_ORD_ID|3\r"), "oids [3]");
        assert_eq!(oids("gnl|BL_ORD_ID|999"), "not found gnl|BL_ORD_ID|999");
        assert_eq!(oids(" gnl|BL_ORD_ID|1"), "not found  gnl|BL_ORD_ID|1");
        assert_eq!(oids("GNL|bl_ord_id|3"), "not found GNL|bl_ord_id|3");
        assert_eq!(oids("3"), "not found 3");
        assert_eq!(oids("myprot9"), "not found myprot9");
        assert_eq!(oids(""), "not found ");
        for tag in ["0", "-1", "abc", "03"] {
            assert_eq!(
                oids(&format!("gnl|BL_ORD_ID|{tag}")),
                "error CObject_id::GetId(): Invalid choice selection: NCBI-General::Object-id.str"
            );
        }
        let parsed = DbIds {
            gi_isam: true,
            str_isam: true,
            ..ids
        };
        assert!(matches!(
            get_oids(&parsed, "gnl|BL_ORD_ID|3"),
            Lookup::Oids(_)
        ));
        assert!(matches!(
            get_oids(&parsed, "3"),
            Lookup::Stop(Stop::Refused(_))
        ));
        assert!(matches!(
            get_oids(&parsed, "XP_000020.1"),
            Lookup::Stop(Stop::Refused(_))
        ));
    }

    #[test]
    fn batch_lines_split_in_two() {
        let dir = std::env::temp_dir().join(format!("losat_bdc_batch_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("b.txt");
        std::fs::write(
            &file,
            "gnl|BL_ORD_ID|3\n  gnl|BL_ORD_ID|4\t3-\t\n\n   \nx 1-5 21\r\n",
        )
        .unwrap();
        let lines = read_batch(file.to_str().unwrap()).unwrap();
        assert_eq!(
            lines,
            [
                ("gnl|BL_ORD_ID|3".to_string(), String::new()),
                ("gnl|BL_ORD_ID|4".to_string(), "3-".to_string()),
                ("x".to_string(), "1-5 21\r".to_string()),
            ]
        );
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(
            matches!(read_batch("in/nosuch.txt"), Err(Stop::Ncbi(message))
            if message == "Argument \"entry_batch\". File is not accessible:  `in/nosuch.txt'")
        );
    }
}
