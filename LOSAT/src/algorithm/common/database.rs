//! What the search programs share in a search of an NCBI BLAST database (`-db NAME`, read
//! through `crate::seqdb`): the opening of the database by its name, the reader's errors
//! as the command line reports them, the `Database:` lines of the pairwise report, and the
//! rules for the stored titles of the reported sequences. BLASTP (stage DB-3), BLASTN
//! (stage DB-4) and TBLASTN (stage DB-5) use them with their own name in LOSAT's refusals;
//! TBLASTN also reads a nucleotide database whole into packed subjects
//! (`read_nucleotide_stores`).
//!
//! The reader (`crate::seqdb`) knows nothing of the command line; this module maps its
//! errors to NCBI's application errors (`crate::cli::NativeError`) and to LOSAT's
//! refusals.

use anyhow::{bail, Result};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};
use crate::sequence::packed_subject::PackedNucleotide;

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:509-535
// ```c++
//     static const string GenerateSearchPath() {
//         string splitter;
//         string path;
// #if defined(NCBI_OS_UNIX)
//         splitter = ":";
// #else
//         splitter = ";";
// #endif
//         // Local directory first;
//         path  = CDirEntry::NormalizePath(CDir::GetCwd(),eFollowLinks);
//         path += splitter;
//         // Then, BLASTDB;
//         CNcbiEnvironment env;
//         path += CDirEntry::NormalizePath(env.Get("BLASTDB"),eFollowLinks);
//         path += splitter;
//         // Finally, the config file.
//         CNcbiApplicationAPI* app = CNcbiApplicationAPI::Instance();
//         if (app) {
//             const CNcbiRegistry& registry = app->GetConfig();
//             if (registry.HasEntry("BLAST", "BLASTDB")) {
//                 path += CDirEntry::NormalizePath(registry.Get("BLAST", "BLASTDB"),eFollowLinks);
//                 path += splitter;
//             }
//         }
//         return path;
//     }
// ```
/// Opens the database `name` (the `-db` value as NCBI keeps it, `seqdb::db_argument_name`)
/// of the sequence type `seqtype` for a search: the one function through which a search
/// program opens a database. The name is looked for as NCBI looks for it: in the current
/// directory, then in the directories of the `BLASTDB` variable, then in those of the
/// registry entry `[BLAST] BLASTDB` (`settings.blastdb`, read from
/// `NCBI_CONFIG__BLAST__BLASTDB`, the program's `.ini` file and `.ncbirc`).
pub(crate) fn open_search_database(
    name: &str,
    seqtype: SeqType,
    settings: &ApplicationSettings,
) -> Result<SeqDb, DatabaseOpenError> {
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => return Err(DatabaseOpenError::Registry(reason.clone())),
    };
    let path = SearchPath::from_environment(registry).map_err(DatabaseOpenError::Reader)?;
    SeqDb::open_in(name, seqtype, &path).map_err(DatabaseOpenError::Reader)
}

/// Why `open_search_database` could not open a database.
#[derive(Debug)]
pub(crate) enum DatabaseOpenError {
    /// The registry entry `[BLAST] BLASTDB` has a value that LOSAT does not read as NCBI
    /// does (`ApplicationSettings::blastdb`).
    Registry(String),
    /// An error of the reader (NCBI's own, or LOSAT's refusal: `database_error`).
    Reader(SeqDbError),
}

/// Opens the database of `-db name` (`open_search_database`) with the command line's
/// errors (`database_error`; LOSAT's refusals name `program`).
pub(crate) fn open_database(
    name: &str,
    seqtype: SeqType,
    settings: &ApplicationSettings,
    program: &str,
) -> Result<SeqDb> {
    open_search_database(name, seqtype, settings).map_err(|error| match error {
        DatabaseOpenError::Registry(reason) => refusal(&reason, program),
        DatabaseOpenError::Reader(error) => database_error(error, program),
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// ```c++
//     catch (const CSeqDBException& e) {                                      \
// 		if (e.GetErrCode() == CSeqDBException::eTooManyOpenFiles) {         \
// 			string err_msg =                                                \
//                 "Too many open files, please raise the open file limit";    \
//         	LOG_POST(Error << "BLAST Database error: " << err_msg);         \
//             exit_code = BLAST_DATABASE_ERROR;                               \
// 		}                                                                   \
// 		else if (e.GetErrCode() == CSeqDBException::eMemoryMappingFailure){ \
// 			string err_msg = "Database memory map file error";              \
//             if (!e.GetMsg().empty())                                        \
//                 err_msg += " - " + e.GetMsg();                              \
//         	LOG_POST(Error << "BLAST Database error: " << err_msg);         \
//             exit_code = BLAST_DATABASE_ERROR;                               \
// 		}                                                                   \
//         else {                                                              \
//         	LOG_POST(Error << "BLAST Database error: " << e.GetMsg());      \
//         	exit_code = BLAST_DATABASE_ERROR;                               \
//         }                                                                   \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:153
// ```c++
// #define BLAST_DATABASE_ERROR                    2
// ```
/// A reader error outside the preliminary search (when the database is opened, and when
/// the sequences and titles of the hit lists are read) as the command line reports it:
/// NCBI's own errors (`SeqDbError::ncbi_message`: a database that is not found, a file
/// that is no database file or cannot be opened) as `BLAST Database error: <text>` with
/// exit status 2; everything else is LOSAT's refusal, in LOSAT's wording (`program` is
/// the program that refuses).
pub(crate) fn database_error(error: SeqDbError, program: &str) -> anyhow::Error {
    match error.ncbi_message() {
        Some(text) => NativeError {
            exit: 2,
            message: format!("BLAST Database error: {text}\n"),
        }
        .into(),
        None => refusal(&error.to_string(), program),
    }
}

/// LOSAT's refusal of a database that `program` does not search as NCBI does.
pub(crate) fn refusal(reason: &str, program: &str) -> anyhow::Error {
    anyhow::anyhow!("{reason}; this database is not supported by LOSAT's {program}")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:240-251
// ```c++
// 	catch (CSeqDBException & e) {
// 		if (e.GetErrCode() == CSeqDBException::eTooManyOpenFiles) {
//             _TRACE(e.GetMsg());
// 			string err_msg = BlastErrorCode2String(BLASTERR_DB_TOO_MANY_OPEN_FILES);
// 			NCBI_THROW(CBlastException, eCoreBlastError, err_msg);
// 		}
// 		else {
//             _TRACE(e.GetMsg());
// 			string err_msg = BlastErrorCode2String(BLASTERR_DB_MEMORY_MAP);
// 			NCBI_THROW(CBlastException, eCoreBlastError, err_msg);
// 		}
// 	}
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// ```c++
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
// ```
/// A reader error while the preliminary search reads the sequences: NCBI's search stage
/// replaces the reader's text (`SeqDbError::ncbi_search_message`), and the application
/// reports `BLAST engine error: <text>` with exit status 3. The error is kept apart from
/// the others (`PreliminaryReadError`), because NCBI has written the outfmt 0 prolog and
/// read the query batch before it.
pub(crate) fn search_read_error(error: SeqDbError, program: &str) -> anyhow::Error {
    match error.ncbi_search_message() {
        Some(text) => PreliminaryReadError(text).into(),
        None => refusal(&error.to_string(), program),
    }
}

/// Bytes of the sequence file that one read of `read_nucleotide_stores` holds at most
/// (`SeqDb::next_block`; a longer sequence is a read of its own). No NCBI counterpart (NCBI
/// maps the volume files); the size does not change what is read.
const STORE_BLOCK_BYTES: u64 = 1 << 20;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1401-1405
// ```c
//     memset((void*) &seq_arg, 0, sizeof(seq_arg));
//
//     /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
//        sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
//     seq_arg.encoding = eBlastEncodingProtein;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:771-775
// ```c
//         } else {
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ```
/// Every sequence of the nucleotide database `db`, in ordinal order, as a packed subject
/// (`PackedNucleotide::from_stored`): the stored 2-bit data, which a translated search's
/// preliminary stage translates as NCBI does (makeblastdb's bases for the ambiguous ones),
/// and the stored ambiguity runs, from which the traceback and the reports read the real
/// codes. The file is read in blocks of at most `STORE_BLOCK_BYTES`, each dropped once its
/// sequences are packed, so the database is held packed only. NCBI reads every sequence's
/// 2-bit data in the first query batch's preliminary search, so a block that cannot be read
/// is that search's error (`search_read_error`).
pub(crate) fn read_nucleotide_stores(db: &SeqDb, program: &str) -> Result<Vec<PackedNucleotide>> {
    let mut stores = Vec::with_capacity(db.num_seqs() as usize);
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, STORE_BLOCK_BYTES) {
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let (length, runs) = stored_parts(&block, i, range.first_ordinal + i as u32, program)?;
            stores.push(PackedNucleotide::from_stored(
                block.nucleotide_packed(i),
                length,
                &runs,
            ));
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(stores)
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the checks of `read_nucleotide_stores` without the copy. The
/// database is read in the same blocks, and every sequence's length and ambiguity runs are read
/// with the same errors in the same order; only the lengths are kept, by ordinal id. A search
/// that reads its subjects where it searches them (`x_read_nucleotide_store`) fails the query
/// batch that the copy would fail, with the same text.
pub(crate) fn x_check_nucleotide_stores(db: &SeqDb, program: &str) -> Result<Vec<usize>> {
    let mut lengths = Vec::with_capacity(db.num_seqs() as usize);
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, STORE_BLOCK_BYTES) {
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let (length, _) = stored_parts(&block, i, range.first_ordinal + i as u32, program)?;
            lengths.push(length);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(lengths)
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): sequence `ordinal` of `db` as `read_nucleotide_stores` packs
/// it (`PackedNucleotide::from_stored` of the stored 2-bit data, length and ambiguity runs), read
/// with one read of the sequence file. The errors are those of `read_nucleotide_stores` for the
/// same sequence. The bytes are read into a buffer that glibc maps from 1 MiB on, taken from
/// `utils::x_mapped::X_BYTES` when it holds one (`SeqDb::x_read_sequence_mapped`), and the store
/// is that buffer cut to the 2-bit data (`PackedNucleotide::x_from_stored_buffer`): a long
/// sequence read for a while and dropped gives its pages back (a store is not given to the pool,
/// which forgets the buffer).
pub(crate) fn x_read_nucleotide_store(
    db: &SeqDb,
    ordinal: u32,
    program: &str,
) -> Result<PackedNucleotide> {
    let block = db
        .x_read_sequence_mapped(ordinal)
        .map_err(|error| search_read_error(error, program))?;
    let (length, runs) = stored_parts(&block, 0, ordinal, program)?;
    let packed = block.x_into_nucleotide_packed();
    if crate::utils::x_mapped::x_is_mapped(packed.capacity()) {
        crate::utils::x_mapped::X_BYTES.x_disown(&packed);
    }
    Ok(PackedNucleotide::x_from_stored_buffer(
        packed, length, &runs,
    ))
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): sequence `i` of `block` (ordinal id `ordinal`), a block of
/// `SeqDb::read_block`, as `read_nucleotide_stores` packs it, with its errors.
pub(crate) fn x_nucleotide_store_in(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    ordinal: u32,
    program: &str,
) -> Result<PackedNucleotide> {
    let (length, runs) = stored_parts(block, i, ordinal, program)?;
    Ok(PackedNucleotide::from_stored(
        block.nucleotide_packed(i),
        length,
        &runs,
    ))
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the length of sequence `i` of `block` (ordinal id
/// `ordinal`), with the errors of `x_nucleotide_store_in`, without the store.
pub(crate) fn x_nucleotide_length_in(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    ordinal: u32,
    program: &str,
) -> Result<usize> {
    stored_parts(block, i, ordinal, program).map(|(length, _)| length)
}

/// The length and the ambiguity runs of sequence `i` of `block` (ordinal id `ordinal`), with the
/// errors of `read_nucleotide_stores` (which reads every sequence through it). NCBI's preliminary
/// search reads the 2-bit data only; it reads the ambiguity data in the traceback, of the
/// reported sequences. LOSAT reads them for every sequence, so data it cannot read are LOSAT's
/// refusal, not an NCBI error.
fn stored_parts(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    ordinal: u32,
    program: &str,
) -> Result<(usize, Vec<crate::seqdb::AmbiguityRun>)> {
    let length = block
        .nucleotide_len(i)
        .map_err(|error| search_read_error(error, program))?;
    let runs = block.ambiguity_runs(i).map_err(|error| {
        refusal(
            &format!("the ambiguity data of database sequence {ordinal} cannot be read ({error})"),
            program,
        )
    })?;
    Ok((length, runs))
}

/// NCBI's error of a database file that cannot be read during the preliminary search
/// (`search_read_error`), with the text after `BLAST engine error: `.
#[derive(Debug)]
pub(crate) struct PreliminaryReadError(pub(crate) &'static str);

impl std::fmt::Display for PreliminaryReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BLAST engine error: {}", self.0)
    }
}

impl std::error::Error for PreliminaryReadError {}

/// A reader error of NCBI's own while the definition lines of the final hit lists are
/// read (`read_title`): `query` is the query whose hit list reached it. NCBI fails that
/// query's batch, after the reports of the batches before.
#[derive(Debug)]
pub(crate) struct HeaderReadError {
    pub(crate) query: usize,
    pub(crate) error: anyhow::Error,
}

impl std::fmt::Display for HeaderReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.error)
    }
}

impl std::error::Error for HeaderReadError {}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:385-391
// ```c++
//     info.name = seqdb->GetDBNameList();
//     info.definition = seqdb->GetTitle();
//     if (info.definition.empty())
//         info.definition = info.name;
//     info.date = seqdb->GetDate();
//     info.total_length = seqdb->GetTotalLength();
//     info.number_seqs = seqdb->GetNumSeqs();
// ```
/// What the `Database:` lines of the pairwise report show of the database `db` opened by
/// the name `name`: the title (the name as given for an empty title), the date and the
/// real totals, also with `-dbsize`.
pub(crate) fn pairwise_database(db: &SeqDb, name: &str) -> PairwiseDatabase {
    let title = db.title();
    PairwiseDatabase {
        definition: if title.is_empty() {
            name.as_bytes().to_vec()
        } else {
            title.to_vec()
        },
        date: db.date().to_string(),
        num_sequences: u64::from(db.num_seqs()),
        total_letters: db.total_length(),
    }
}

/// The stored title of a database sequence, as the reports read it (`read_title`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct StoredTitle {
    /// The title of the sequence's definition line, with a byte outside ' ' ..= '~' read
    /// as '#'.
    pub(crate) title: Vec<u8>,
    /// Whether the stored title held a byte outside ' ' ..= '~' (`check_reported_title`).
    pub(crate) replaced_bytes: bool,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310
// ```c++
// CTempString
// CSeqDBVol::x_GetHdrAsn1Binary(int oid) const
// {
//     TIndx hdr_start = 0;
//     TIndx hdr_end   = 0;
//
//     //m_Atlas.Lock(locked);
//
//     if (!m_HdrFileOpened) x_OpenHdrFile();
// ```
// NCBI reference for the titles:
// c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:970-1017 (`s_GetBioseqTitle`: the
// title of the first definition line, the others appended as ` >id title`).
/// Reads the definition line of the database sequence `ordinal` for a report, where the
/// hit list of query `query` reaches it. A reader error of NCBI's own is a
/// `HeaderReadError` that names the query (NCBI fails that query's batch); LOSAT's
/// refusals (`program` names the program) are errors of their own.
pub(crate) fn read_title(
    db: &SeqDb,
    ordinal: u32,
    query: usize,
    program: &str,
) -> Result<StoredTitle> {
    stored_title(db.deflines(ordinal), ordinal, query, program)
}

/// `read_title` of the sequence `ordinal` whose definition lines `deflines` were read
/// (`SeqDb::deflines`, or `SeqDb::visit_deflines` with the same results).
pub(crate) fn stored_title(
    deflines: std::result::Result<Vec<crate::seqdb::OrdinalDefline>, crate::seqdb::SeqDbError>,
    ordinal: u32,
    query: usize,
    program: &str,
) -> Result<StoredTitle> {
    let mut deflines = deflines.map_err(|error| {
        if error.ncbi_message().is_some() {
            anyhow::Error::new(HeaderReadError {
                query,
                error: database_error(error, program),
            })
        } else {
            database_error(error, program)
        }
    })?;
    // makeblastdb without -parse_seqids writes one definition line per sequence; NCBI
    // shows the others of a sequence as ` >id title` after the first.
    if deflines.len() != 1 {
        bail!(
            "database sequence {ordinal} has {} definition lines; a sequence with more or fewer than one is not supported by LOSAT's {program}",
            deflines.len()
        );
    }
    let defline = deflines.remove(0);
    Ok(StoredTitle {
        title: defline.title,
        replaced_bytes: defline.replaced_title_bytes > 0,
    })
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1865-1886
// ```c++
//         string message;
//         if (io != NULL) {
//             message += io->GetStackTrace() + "\n";
//         }
//         message += "Bad char [0x" +
//                    NStr::NumericToString((unsigned char)c,0,16)+
//                    "] in string";
//         if (io != NULL) {
//             message += " at " + io->GetPosition();
//         }
//         if (!str.empty()) {
//             message += "\n" + str;
//         }
//         switch (fix_method) {
//         case eFNP_ReplaceAndWarn:
// #if 0
//             CNcbiDiag(eDiag_Error, eDPF_Default)
//                 << ErrCode(NCBI_ERRCODE_X, 7) << message << Endm;
// #else
//             ERR_POST_X(7, message);
// #endif
//             break;
// ```
// NCBI BLAST+ 2.17.0 `blastp -db` on a version 4 database whose header file holds
// the bytes 0xE9 and 0x7F in two titles (campaign dbsearch-20261012,
// `$BUILD_ROOT/dbsearch/db3/exp_badbyte`): the report shows `#` for the byte, and
// standard error gets, for each reported subject and bad byte,
//   Error: (803.7) [blastp] Blast-def-line-set.E.title
//   Bad char [0xE9] in string at byte 21
//   <the stored title>
// as often as the formatter reads the definition line: twice per subject for the
// default outfmt 6, three times for `6 sseqid stitle`, five times for outfmt 0. The
// number of readings is the formatter's, not ported; makeblastdb and `losat makedb`
// store `#` for such a byte, so only a header file written otherwise holds one.
/// LOSAT's rejection of the reported database sequence `ordinal` (one with rows in the
/// final hit lists) whose stored title held a byte outside ' ' ..= '~' (`replaced_bytes`,
/// `StoredTitle`): NCBI writes an error message for the byte each time its formatter
/// reads the definition line.
pub(crate) fn check_reported_title(
    replaced_bytes: bool,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    if replaced_bytes {
        bail!(
            "database sequence {ordinal} has a stored title with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line); this is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of the database
/// sequence `ordinal` with the stored title `title` (`report::defline::database_subject_id`;
/// `protein` for a protein database). A title whose id LOSAT cannot write as NCBI is
/// LOSAT's refusal (`program` names the program).
pub(crate) fn tabular_subject_id(
    title: &[u8],
    ordinal: u32,
    protein: bool,
    program: &str,
) -> Result<Vec<u8>> {
    crate::report::defline::database_subject_id(title, protein).map_err(|problem| {
        anyhow::anyhow!(
            "database sequence {ordinal} has a title {problem}, which is not supported by LOSAT's {program}"
        )
    })
}

/// LOSAT's rejection of the outfmt 0 title `title` of the protein database sequence
/// `ordinal` that it does not write as NCBI: a title of punctuation that NCBI's
/// `x_CleanAndCompress` reads past the end of when it makes the title of the description
/// table (NCBI BLAST+ 2.17.0 crashes there), as for the records of a subject set
/// (`check_shown_subject_titles`). Approved exception 2 of PD-LOSAT-NCBI-DEFECTS covers
/// the nucleotide subjects of BLASTN, so BLASTN does not use it.
pub(crate) fn check_shown_protein_title(title: &[u8], ordinal: u32, program: &str) -> Result<()> {
    let reads_past_end = crate::report::defline::generate_defline(title, true, true)
        .is_ok_and(|defline| defline.reads_past_end);
    if reads_past_end {
        bail!(
            "database sequence {ordinal} has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0); this is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// Byte 0 of NCBIstdaa: the gap character (`-`), and the sentinel byte between the
/// sequences of a protein volume's sequence file.
const NCBISTDAA_GAP: u8 = 0;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593
// ```c
//             /* skip special characters */
//             if (i == AMINOACID_TO_NCBISTDAA['U'] ||
//                 i == AMINOACID_TO_NCBISTDAA['O'] ||
//                 i == AMINOACID_TO_NCBISTDAA['-'] ||
// ```
// A stored gap character (0, which makeblastdb writes for `-` in a protein FASTA)
// never aligns in NCBI (row and column 0 stay `BLAST_SCORE_MIN`) and NCBI's formatter
// then counts it as a gap of the alignment; LOSAT's kernels read byte 0 as the end of a
// sequence, so a sequence that holds one is rejected. A byte that is no NCBIstdaa
// residue (28 or more) indexes NCBI's matrix past its rows and is rejected too.
/// LOSAT's refusal of the protein database sequence `ordinal` (`with_sentinels`: its
/// residues with the byte before and the byte after them in the sequence file) when
/// `program` cannot search it as NCBI does: a sequence not between two 0 bytes, one
/// that holds the gap character, or one with a byte that is no NCBIstdaa residue.
pub(crate) fn check_protein_sequence(
    with_sentinels: &[u8],
    ordinal: u64,
    program: &str,
) -> Result<()> {
    let (first, last) = (with_sentinels[0], with_sentinels[with_sentinels.len() - 1]);
    if first != NCBISTDAA_GAP || last != NCBISTDAA_GAP {
        bail!(
            "database sequence {ordinal} is not between two NUL bytes in its sequence file (the volume is damaged), which is not supported by LOSAT's {program}"
        );
    }
    let residues = &with_sentinels[1..with_sentinels.len() - 1];
    // Every residue must be a byte other than the gap and below BLASTAA_SIZE: with the gap byte
    // 0 (asserted below), the residues 1..BLASTAA_SIZE, the bytes whose `residue - 1` (wrapping,
    // so the gap becomes 255) is below BLASTAA_SIZE - 1. The largest such difference, taken over
    // all bytes without an early exit, is a loop the compiler vectorizes; the early exit of
    // `all` made it a loop over single bytes (LOSAT-only cost; NCBI checks nothing here).
    const _: () = assert!(
        NCBISTDAA_GAP == 0
            && crate::utils::matrix::BLASTAA_SIZE > 1
            && crate::utils::matrix::BLASTAA_SIZE <= 256
    );
    const FIRST: u8 = NCBISTDAA_GAP + 1;
    const LETTERS: usize = crate::utils::matrix::BLASTAA_SIZE - FIRST as usize;
    let largest = residues
        .iter()
        .map(|&residue| residue.wrapping_sub(FIRST))
        .max()
        .unwrap_or(0);
    if usize::from(largest) < LETTERS {
        return Ok(());
    }
    if let Some(position) = residues
        .iter()
        .position(|&residue| residue == NCBISTDAA_GAP)
    {
        bail!(
            "database sequence {ordinal} holds a gap character ('-') at position {}; a protein database sequence with gap characters is not supported by LOSAT's {program} (NCBI BLAST+ never aligns the character and then counts it as a gap of the alignment)",
            position + 1
        );
    }
    if let Some(position) = residues
        .iter()
        .position(|&residue| usize::from(residue) >= crate::utils::matrix::BLASTAA_SIZE)
    {
        bail!(
            "database sequence {ordinal} holds a byte that is no NCBIstdaa residue at position {}, which is not supported by LOSAT's {program}",
            position + 1
        );
    }
    Ok(())
}

/// Bytes of the sequence file that one read of `read_protein_stores` holds at most
/// (`SeqDb::next_block`). No NCBI counterpart: NCBI maps the volume files. The size
/// changes nothing but the number of reads.
const PROTEIN_STORE_BLOCK_BYTES: u64 = 1 << 20;

/// The sequences of a protein database as its sequence files store them: NCBIstdaa
/// bytes, each sequence followed by the sentinel byte 0 (and the first preceded by
/// one), in ordinal order, with the totals of the index files. NCBI searches these
/// bytes in place (blast_engine.c:1404, `eBlastEncodingProtein`), and its traceback and
/// composition redo read them again; a program that holds them all searches, redoes and
/// shows the bytes NCBI does.
#[derive(Debug)]
pub(crate) struct StoredProteins {
    /// The sentinel, then each sequence and the sentinel after it.
    bytes: Vec<u8>,
    /// Where each sequence's residues start in `bytes`, and one more entry: sequence
    /// `i` is `bytes[starts[i]..starts[i + 1] - 1]`.
    starts: Vec<usize>,
    /// The total length of the index files (`CSeqDB::GetTotalLength`).
    total_length: u64,
}

impl StoredProteins {
    /// The number of sequences.
    pub(crate) fn len(&self) -> usize {
        self.starts.len() - 1
    }

    /// Whether the database holds no sequence.
    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The length of sequence `ordinal`.
    pub(crate) fn length(&self, ordinal: usize) -> usize {
        self.starts[ordinal + 1] - self.starts[ordinal] - 1
    }

    /// The residues of sequence `ordinal` followed by the sentinel byte 0.
    pub(crate) fn residues_with_sentinel(&self, ordinal: usize) -> &[u8] {
        &self.bytes[self.starts[ordinal]..self.starts[ordinal + 1]]
    }

    /// EXPERIMENT (LOSAT_X_BXUNGAPSENT): `residues_with_sentinel` with the sentinel before the
    /// residues (the store's first byte, or the sentinel after the sequence before).
    pub(crate) fn x_residues_with_sentinels(&self, ordinal: usize) -> &[u8] {
        &self.bytes[self.starts[ordinal] - 1..self.starts[ordinal + 1]]
    }

    /// The total length of the index files.
    pub(crate) fn total_length(&self) -> u64 {
        self.total_length
    }

    /// The length of the longest sequence (0 for an empty database).
    pub(crate) fn max_length(&self) -> usize {
        (0..self.len()).map(|i| self.length(i)).max().unwrap_or(0)
    }
}

/// Reads every sequence of the protein database `db` once, in blocks of
/// `PROTEIN_STORE_BLOCK_BYTES` in ordinal order, as its sequence files store them. A
/// reader error is NCBI's error of the preliminary search (`search_read_error`: NCBI
/// reads every sequence there); a sequence `program` cannot search as NCBI does is
/// LOSAT's refusal (`check_protein_sequence`).
pub(crate) fn read_protein_stores(db: &SeqDb, program: &str) -> Result<StoredProteins> {
    let sequences = db.num_seqs() as usize;
    let mut bytes = Vec::with_capacity(
        usize::try_from(db.total_length())
            .unwrap_or(0)
            .saturating_add(sequences + 1),
    );
    let mut starts = Vec::with_capacity(sequences + 1);
    bytes.push(NCBISTDAA_GAP);
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, PROTEIN_STORE_BLOCK_BYTES) {
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let with_sentinels = block.protein_with_sentinels(i);
            check_protein_sequence(
                with_sentinels,
                u64::from(range.first_ordinal) + i as u64,
                program,
            )?;
            starts.push(bytes.len());
            bytes.extend_from_slice(&with_sentinels[1..]);
        }
        next_ordinal = range.end_ordinal();
    }
    starts.push(bytes.len());
    Ok(StoredProteins {
        bytes,
        starts,
        total_length: db.total_length(),
    })
}

/// EXPERIMENT (LOSAT_X_BXDBSTREAM): `read_protein_stores` without the store: the same reads in
/// the same blocks and the same checks in the same order, with the same errors, the number of
/// residues of each sequence handed to `each` (ordinal id, residues) and its bytes dropped with
/// its block. A search that reads the sequences again later (BLASTX's passes,
/// `blastx/x_db_stream.rs`) knows from it that they read, and what the store's
/// `max_length` would be.
pub(crate) fn x_check_protein_stores(
    db: &SeqDb,
    program: &str,
    mut each: impl FnMut(u32, usize),
) -> Result<()> {
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, PROTEIN_STORE_BLOCK_BYTES) {
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let with_sentinels = block.protein_with_sentinels(i);
            let ordinal = range.first_ordinal + i as u32;
            check_protein_sequence(with_sentinels, u64::from(ordinal), program)?;
            each(ordinal, with_sentinels.len() - 2);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::matrix::BLASTAA_SIZE;

    /// The residues between the two sentinel bytes.
    fn stored(residues: &[u8]) -> Vec<u8> {
        let mut bytes = vec![NCBISTDAA_GAP];
        bytes.extend_from_slice(residues);
        bytes.push(NCBISTDAA_GAP);
        bytes
    }

    fn check(residues: &[u8]) -> std::result::Result<(), String> {
        check_protein_sequence(&stored(residues), 7, "blastx").map_err(|e| e.to_string())
    }

    /// The check before the single pass (an `all` over the residues).
    fn accepted_before(residues: &[u8]) -> bool {
        residues
            .iter()
            .all(|&residue| residue != NCBISTDAA_GAP && usize::from(residue) < BLASTAA_SIZE)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593 (quoted at
    // `check_protein_sequence`): the gap byte and bytes past the matrix are what the check
    // refuses. Not a port of an NCBI unit case.
    #[test]
    fn the_check_accepts_the_letters_and_names_the_first_bad_byte() {
        let first = NCBISTDAA_GAP + 1;
        let last = (BLASTAA_SIZE - 1) as u8;
        assert_eq!(check(&[]), Ok(()));
        assert_eq!(check(&[first]), Ok(()));
        assert_eq!(check(&[last]), Ok(()));
        let letters: Vec<u8> = (0..1_000).map(|i| [first, last][i % 2]).collect();
        assert_eq!(check(&letters), Ok(()));
        // (byte, words of the refusal, the character after its position)
        for (bad, says, after) in [
            (NCBISTDAA_GAP, "gap character", ';'),
            (BLASTAA_SIZE as u8, "no NCBIstdaa residue", ','),
            (255, "no NCBIstdaa residue", ','),
        ] {
            for length in [1, 2, 31, 32, 33, 100, 1_000] {
                for at in [0, length / 2, length - 1] {
                    let mut residues = letters[..length].to_vec();
                    residues[at] = bad;
                    let error = check(&residues).unwrap_err();
                    assert!(
                        error.contains(says)
                            && error.contains(&format!("position {}{after}", at + 1)),
                        "byte {bad} at {at} of {length}: {error}"
                    );
                }
            }
        }
        // A gap before a byte past the matrix is reported as the gap, as before.
        let mut residues = letters[..10].to_vec();
        residues[7] = NCBISTDAA_GAP;
        residues[3] = 200;
        assert!(check(&residues).unwrap_err().contains("gap character"));
        // The sentinels.
        assert!(check_protein_sequence(&[1, first, 0], 7, "blastx").is_err());
        assert!(check_protein_sequence(&[0, first, 1], 7, "blastx").is_err());
    }

    // Random residues around the letters' bounds: the single pass accepts exactly what the
    // `all` accepted.
    #[test]
    fn the_check_accepts_what_the_check_before_it_accepted() {
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for case in 0..2_000 {
            let length = (next() % 200) as usize;
            let mut residues: Vec<u8> = (0..length).map(|_| 1 + (next() % 27) as u8).collect();
            if case % 2 == 0 && length > 0 {
                let at = (next() as usize) % length;
                residues[at] = [0, 27, 28, 29, 254, 255][(next() % 6) as usize];
            }
            assert_eq!(
                check(&residues).is_ok(),
                accepted_before(&residues),
                "{residues:?}"
            );
        }
    }
}

/// The taxonomy database of a run (`CTaxDBFileInfo`, a static of NCBI's `GetTaxNames`) and the
/// taxids the volumes have looked up (`CSeqDBVol::x_GetTaxonomy` caches them per volume).
///
/// It lives on the database of the search, so a `run_local` call that writes outfmt 15 and 16
/// of one `-db` search shares it: the open and the volume-lookup messages are written for the
/// first of the two only (the CLI writes one format). NCBI, which runs each format alone, writes
/// them for both; the ABI callers that ask for both do not see them twice.
#[derive(Default)]
pub(crate) struct TaxdbRun {
    lazy: std::sync::Mutex<Option<std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>>>,
    volume_lookups: std::sync::Mutex<std::collections::HashSet<(usize, i32)>>,
}

impl TaxdbRun {
    /// The run's taxonomy database, a handle that opens it at its first name asked for.
    fn taxdb(
        &self,
        search_path: &crate::seqdb::SearchPath,
    ) -> std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb> {
        let mut lazy = self.lazy.lock().expect("the taxonomy database");
        std::sync::Arc::clone(lazy.get_or_insert_with(|| {
            std::sync::Arc::new(crate::seqdb::taxdb::LazyTaxDb::new(search_path.clone()))
        }))
    }

    /// EXPERIMENT (LOSAT_X_BPQBLOCKFMTSHADOW): forgets the opened taxonomy database and the taxids
    /// looked up, so that a second search of the same queries (the shadow's search of all queries)
    /// writes the messages of opening and of the volumes' lookups as the first search did.
    pub(crate) fn x_reset(&self) {
        *self.lazy.lock().expect("the taxonomy database") = None;
        self.volume_lookups
            .lock()
            .expect("the volume lookups")
            .clear();
    }
}

/// LOSAT's refusal of a taxonomy database name or record (`read_reported_taxids`): a scientific
/// name NCBI's XML stream cannot write, or a record NCBI reads from unmapped memory. A run that
/// holds its XML2 or JSON document drops it on this error and no other (BLASTX,
/// `XmlRefusal::of_taxonomy`); the other errors of the lookups (a read error of a file) are not
/// refusals, and the held document is written as NCBI writes what it has.
#[derive(Debug)]
pub(crate) struct TaxonomyRefusal(String);

impl std::fmt::Display for TaxonomyRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TaxonomyRefusal {}

/// [`taxdb_error`] for the tests of the callers' error handling.
#[cfg(test)]
pub(crate) fn taxdb_error_for_tests(error: SeqDbError) -> anyhow::Error {
    taxdb_error(error, "BLASTX")
}

/// The error of a taxonomy database lookup or open: [`TaxonomyRefusal`] for a refusal of the
/// reader (`SeqDbError::Unsupported`), [`database_error`] for the others.
fn taxdb_error(error: SeqDbError, program: &str) -> anyhow::Error {
    let refused = matches!(error, SeqDbError::Unsupported { .. });
    let error = database_error(error, program);
    if refused {
        anyhow::Error::new(TaxonomyRefusal(format!("{error:#}")))
    } else {
        error
    }
}

/// The Blast-def-line taxid of each reported subject by rank (0 for none and for a subject
/// that is not in `reported`), which the XML2 and JSON formats write in the `Hit` description,
/// and what the taxonomy database says about it. They are read apart from the titles, for those
/// formats only. `ordinals[rank]` is the ordinal id of the subject of rank `rank`.
///
/// NCBI writes the scientific name of a taxid when it finds `taxdb.bti` on the search path
/// (`data4xml2format.cpp x_InitCommon`, `blastxml2_format.cpp:380-386`; validity is not
/// checked): `taxonomy[rank]` is then `Named`, or `NotFound` when the taxonomy database does
/// not know the taxid (the description is the id alone). Without `taxdb.bti`, or without a
/// taxid, `taxonomy` is empty. A scientific name that NCBI's XML stream cannot write, and a
/// taxonomy database that LOSAT cannot read as NCBI does, are refusals made here, before any
/// byte is written.
///
/// NCBI's messages about the taxonomy database go to standard error, and not here: those of
/// opening and of the lookup that the volume makes for each taxid of a fetched subject
/// (`CSeqDBVol::x_GetTaxonomy`, once per volume and taxid) are in the returned [`TaxdbNotes`],
/// which writes them at the first `Hit` that needs them, and those of the lookup of each `Hit`
/// are in `taxonomy` and written at the hit ([`TaxdbPrinter::at_hit`]). `program` is the
/// program that refuses, as LOSAT's refusals name it (`BLASTP`); its messages name it in lower
/// case.
pub(crate) fn read_reported_taxids(
    db: &SeqDb,
    run: &TaxdbRun,
    ordinals: &[u32],
    reported: &std::collections::BTreeSet<usize>,
    program: &str,
) -> Result<(
    Vec<i32>,
    Vec<crate::report::blast_xml::HitTaxonomy>,
    TaxdbNotes,
)> {
    use crate::report::blast_xml::HitTaxonomy;
    let mut taxids = vec![0; ordinals.len()];
    for &rank in reported {
        let defline_taxids = db
            .defline_taxids(ordinals[rank])
            .map_err(|error| database_error(error, program))?;
        // One definition line per reported sequence (`stored_title`).
        taxids[rank] = defline_taxids.first().copied().unwrap_or(0);
    }
    if taxids.iter().all(|&taxid| taxid == 0) || !db.taxdb_found() {
        return Ok((taxids, Vec::new(), TaxdbNotes::default()));
    }
    let Some(search_path) = db.search_path() else {
        return Ok((taxids, Vec::new(), TaxdbNotes::default()));
    };
    let taxdb = run.taxdb(search_path);
    let mut notes = TaxdbNotes {
        program: program.to_ascii_lowercase(),
        ..TaxdbNotes::default()
    };
    // NCBI opens the database at its first lookup; the texts are written at the first `Hit`.
    notes.open = taxdb.open().map_err(|error| taxdb_error(error, program))?;
    let mut asked: std::collections::BTreeMap<i32, HitTaxonomy> = Default::default();
    let mut taxonomy = vec![HitTaxonomy::Unnamed; taxids.len()];
    notes.volumes = vec![0; taxids.len()];
    for &rank in reported {
        let taxid = taxids[rank];
        if taxid == 0 {
            continue;
        }
        // The volume's lookup when it fetches the sequence (once per volume and taxid;
        // `x_GetTaxonomy` asks for a positive taxid only, seqdbvol.cpp:1136).
        let volume = db
            .volume_index(ordinals[rank])
            .map_err(|error| database_error(error, program))?;
        notes.volumes[rank] = volume;
        let first = taxid > 0
            && run
                .volume_lookups
                .lock()
                .expect("the volume lookups")
                .insert((volume, taxid));
        if first {
            let mut messages = Vec::new();
            taxdb
                .names(taxid, &mut messages)
                .map_err(|error| taxdb_error(error, program))?;
            notes.lookups.insert((volume, taxid), messages);
        }
        // The lookup of the `Hit` (`CSeqDB::GetTaxInfo`).
        if !asked.contains_key(&taxid) {
            let mut messages = Vec::new();
            let found = taxdb
                .names(taxid, &mut messages)
                .map_err(|error| taxdb_error(error, program))?;
            let answer = match found {
                Some(names) => {
                    crate::report::blast_xml::names::check_text(&names.scientific_name).map_err(
                        |refusal| {
                            anyhow::Error::new(TaxonomyRefusal(
                                refusal
                                    .error(program, "a taxonomy database scientific name")
                                    .to_string(),
                            ))
                        },
                    )?;
                    HitTaxonomy::Named(names.scientific_name)
                }
                None => HitTaxonomy::NotFound(messages),
            };
            asked.insert(taxid, answer);
        }
        taxonomy[rank] = asked[&taxid].clone();
    }
    Ok((taxids, taxonomy, notes))
}

/// How many taxonomy database entries cost as much to check as one database sequence: the
/// taxonomy side is chosen when `entries <= TAXDB_ENTRIES_PER_SEQUENCE * sequences`. Measured
/// (OF5b fix, 1 thread, hot page cache): the pass over NCBI's taxdb takes 0.15 s for 2,314,695
/// entries (about 0.065 us each), the database side 0.25 s for 2,000,000 sequences (about
/// 0.125 us each), so an entry costs about half a sequence.
const TAXDB_ENTRIES_PER_SEQUENCE: u64 = 2;

/// The most database sequences the per-sequence side of [`taxdb_may_refuse`] reads after the
/// taxonomy side found a record that could be refused (about 0.3 us each, so a few seconds);
/// a larger database is held without reading it.
const PREPASS_DB_SEQUENCES_MAX: u64 = 1 << 22;

/// Whether [`read_reported_taxids`] could refuse for some sequence of the database, which a
/// run that writes its batches one by one cannot take back (OF5b BLASTX, DECISIONS 32): a
/// scientific name NCBI's XML stream cannot write, a record past the last mapped page of
/// `taxdb.btd`, or a header or the taxonomy database that cannot be read. A superset of what a
/// search reports: limits (`-gilist`, `-taxidlist`) and the hits themselves are ignored.
/// `false` without `taxdb.bti` on the search path and for a database without taxids.
///
/// The answer comes from whichever side is cheaper: the taxonomy side when the entries of
/// `taxdb.bti` (from its header) number at most [`TAXDB_ENTRIES_PER_SEQUENCE`] times the
/// sequences of the database (`SeqDb::num_seqs`), else the database side:
/// * [`taxdb_side`]: one sequential pass over the taxonomy database asks whether any taxid's
///   lookup could be refused; none -> `false`, without reading the database. This is the
///   normal case for NCBI's taxdb, so the cost does not grow with the database.
/// * [`database_side`]: the first definition line's taxid of every sequence, each distinct
///   non-zero one looked up. It decides the cases the taxonomy side leaves open. When the
///   taxonomy side ran first (the taxonomy database has at most twice as many entries as the
///   database has sequences) and found something, a database of more than
///   [`PREPASS_DB_SEQUENCES_MAX`] sequences is held without this pass.
///
/// "No hold" means that no lookup of the run can refuse after a byte: the taxonomy side proves it
/// for every taxid, the database side for every taxid the database holds.
///
/// It has its own taxonomy handle (`LazyTaxDb`), so the run's [`TaxdbRun`], its volume lookups
/// and the notes that [`TaxdbPrinter`] writes are as they would be without this call. That
/// handle and the run's are two opens of the same files, which differ only if `taxdb.bti` or
/// `taxdb.btd` is replaced while the run goes on (the run opens its handle at its first lookup,
/// possibly many batches after this call): the run can then refuse after bytes were written.
/// NCBI reads the two files at its first lookup too, and a replaced file gives it a mixed
/// result as well; this is not guarded (DECISIONS 32).
pub(crate) fn taxdb_may_refuse(db: &SeqDb) -> bool {
    use crate::seqdb::taxdb::LazyTaxDb;
    if !db.taxdb_found() {
        return false;
    }
    let Some(search_path) = db.search_path() else {
        return false;
    };
    let taxdb = LazyTaxDb::new(search_path.clone());
    let sequences = u64::from(db.num_seqs());
    // An index that cannot be opened is left to the database side, which asks for a name only
    // when a sequence has a taxid.
    let taxdb_first = taxdb.entry_count().map_or(false, |entries| {
        entries <= TAXDB_ENTRIES_PER_SEQUENCE * sequences
    });
    if taxdb_first {
        if taxdb_side(&taxdb) == Some(false) {
            return false;
        }
        if sequences > PREPASS_DB_SEQUENCES_MAX {
            return true;
        }
    }
    database_side(db, &taxdb)
}

/// The taxonomy side of [`taxdb_may_refuse`]: `Some(false)` when no lookup of any taxid can be
/// refused or give a scientific name that `names::check_text` refuses
/// ([`crate::seqdb::taxdb::TaxDb::any_lookup_may_refuse`]); `Some(true)` when one can; `None` when
/// the taxonomy database cannot be opened or read, which only the lookups of the database's taxids
/// can answer.
fn taxdb_side(taxdb: &crate::seqdb::taxdb::LazyTaxDb) -> Option<bool> {
    taxdb
        .any_lookup_may_refuse(&|name| crate::report::blast_xml::names::check_text(name).is_err())
        .ok()
}

/// The database side of [`taxdb_may_refuse`]: whether the lookup of the first definition line's
/// taxid of some sequence would be refused, looking each distinct non-zero taxid up with the
/// `names` + `names::check_text` code `read_reported_taxids` uses. The headers are read in
/// ranges of [`crate::seqdb::HEADER_RANGE_SEQS`] sequences, one read per MiB
/// ([`SeqDb::visit_defline_taxids`]). A header that cannot be read answers `true`. Memory: the
/// distinct taxids, nothing per sequence.
fn database_side(db: &SeqDb, taxdb: &crate::seqdb::taxdb::LazyTaxDb) -> bool {
    let mut seen: std::collections::HashSet<i32> = Default::default();
    let mut refused = false;
    for range in db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS) {
        let whole = db.visit_defline_taxids(range, |_ordinal, defline_taxids| {
            let Ok(defline_taxids) = defline_taxids else {
                refused = true;
                return false;
            };
            // One definition line per sequence (`read_reported_taxids`).
            let taxid = defline_taxids.first().copied().unwrap_or(0);
            if taxid == 0 || !seen.insert(taxid) {
                return true;
            }
            let mut messages = Vec::new();
            refused = match taxdb.names(taxid, &mut messages) {
                Err(_) => true,
                Ok(Some(names)) => {
                    crate::report::blast_xml::names::check_text(&names.scientific_name).is_err()
                }
                Ok(None) => false,
            };
            !refused
        });
        if !whole {
            return refused;
        }
    }
    false
}

/// NCBI's messages about the taxonomy database that are not the lookup of a `Hit`: the texts of
/// opening it (at the first lookup of the run) and the lookup that a volume makes the first
/// time it fetches a subject with a taxid (`CSeqDBVol::x_GetTaxonomy`). NCBI writes them when
/// the first `Hit` that needs them is formatted, after the reader's messages of its batch, so
/// they are kept here and written by [`TaxdbPrinter::at_hit`].
#[derive(Default)]
pub(crate) struct TaxdbNotes {
    /// The program's name in lower case, as the messages show it.
    program: String,
    /// The texts of opening the database, when this batch opened it.
    open: Vec<String>,
    /// The volume of each rank with a taxid.
    volumes: Vec<usize>,
    /// The texts of the lookup of each (volume, taxid) that no earlier batch had looked up.
    lookups: std::collections::BTreeMap<(usize, i32), Vec<String>>,
}

impl TaxdbNotes {
    /// A printer that writes these notes as the `Hit`s of a report come.
    pub(crate) fn printer(&self) -> TaxdbPrinter<'_> {
        TaxdbPrinter {
            notes: self,
            opened: false,
            written: Default::default(),
        }
    }
}

/// Writes the taxonomy database's messages of a report in NCBI's order: the notes of the
/// taxid of a `Hit` before the messages of the `Hit`'s own lookup.
pub(crate) struct TaxdbPrinter<'a> {
    notes: &'a TaxdbNotes,
    opened: bool,
    written: std::collections::HashSet<(usize, i32)>,
}

impl TaxdbPrinter<'_> {
    /// NCBI's messages at the `Hit` of subject `rank` (taxid `taxid`, `taxonomy` what its lookup
    /// said): the texts of opening and of the volume's lookup, where they first apply, then the
    /// messages of the `Hit`'s lookup (`CSeqDB::GetTaxInfo`, every time), to `sink`.
    pub(crate) fn at_hit(
        &mut self,
        rank: usize,
        taxid: i32,
        taxonomy: &crate::report::blast_xml::HitTaxonomy,
        sink: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        use crate::seqdb::taxdb::err_post_line;
        let mut text = String::new();
        let program = &self.notes.program;
        let mut push = |line: &str| {
            text.push_str(&err_post_line(program, line));
            text.push('\n');
        };
        if taxid != 0 {
            if !std::mem::replace(&mut self.opened, true) {
                self.notes.open.iter().for_each(|line| push(line));
            }
            if let Some(volume) = self.notes.volumes.get(rank).copied() {
                if let Some(lines) = self.notes.lookups.get(&(volume, taxid)) {
                    if self.written.insert((volume, taxid)) {
                        lines.iter().for_each(|line| push(line));
                    }
                }
            }
        }
        if let crate::report::blast_xml::HitTaxonomy::NotFound(messages) = taxonomy {
            messages.iter().for_each(|line| push(line));
        }
        if text.is_empty() {
            return Ok(());
        }
        sink.write_all(text.as_bytes())?;
        sink.flush()
    }
}

#[cfg(test)]
mod taxdb_notes_tests {
    use super::*;
    use crate::report::blast_xml::HitTaxonomy;

    fn lines(printer: &mut TaxdbPrinter<'_>, rank: usize, taxid: i32, hit: &HitTaxonomy) -> String {
        let mut sink = Vec::new();
        printer.at_hit(rank, taxid, hit, &mut sink).unwrap();
        String::from_utf8(sink).unwrap()
    }

    /// NCBI writes the texts of opening at the first `Hit` with a taxid, the volume's lookup
    /// at the first `Hit` of each (volume, taxid), and the `Hit`'s own lookup every time, each
    /// after the lines before it (and so after the reader's messages of the batch).
    #[test]
    fn the_notes_come_at_the_first_hit_that_needs_them() {
        let notes = TaxdbNotes {
            program: "blastp".to_string(),
            open: vec!["Error: open".to_string()],
            volumes: vec![0, 0, 1, 0],
            lookups: [
                ((0, 9606), vec!["Error: lookup 9606".to_string()]),
                ((1, 9606), vec!["Error: lookup 9606".to_string()]),
            ]
            .into_iter()
            .collect(),
        };
        let mut printer = notes.printer();
        let missing = HitTaxonomy::NotFound(vec!["Error: hit".to_string()]);
        // a subject without a taxid writes nothing and does not open the database
        assert_eq!(lines(&mut printer, 3, 0, &HitTaxonomy::Unnamed), "");
        assert_eq!(
            lines(&mut printer, 0, 9606, &missing),
            "Error: [blastp] Error: open\nError: [blastp] Error: lookup 9606\nError: [blastp] Error: hit\n"
        );
        // the same (volume, taxid): the volume does not look it up again, the `Hit` does
        assert_eq!(
            lines(&mut printer, 1, 9606, &missing),
            "Error: [blastp] Error: hit\n"
        );
        assert_eq!(
            lines(&mut printer, 1, 9606, &HitTaxonomy::Named(Vec::new())),
            ""
        );
        // another volume: its own lookup
        assert_eq!(
            lines(&mut printer, 2, 9606, &HitTaxonomy::Named(Vec::new())),
            "Error: [blastp] Error: lookup 9606\n"
        );
        // no notes (a batch after the first, or a database without a taxonomy database)
        let none = TaxdbNotes::default();
        assert_eq!(
            lines(&mut none.printer(), 0, 9606, &HitTaxonomy::Unnamed),
            ""
        );
    }
}

#[cfg(test)]
mod taxdb_prepass_tests {
    //! `taxdb_may_refuse`: the answer of the lookup of every sequence's taxid, with the
    //! databases of `tests/fixtures/outfmt16/db` (`sprot_tax`: two sequences, taxid 9606;
    //! `sprot_mixtax`: those two with taxid 10090 and a third with taxid 9606) and the taxonomy
    //! databases of `tests/seqdb/taxdb`.
    use super::*;
    use crate::seqdb::alias::SearchPath;
    use crate::seqdb::SeqType;
    use std::path::{Path, PathBuf};

    fn root() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")
    }

    /// The database `name` of `tests/fixtures/outfmt16/db` opened with `BLASTDB=taxdb`
    /// (`None`: not set), from an empty current directory.
    fn database(name: &str, taxdb: Option<&Path>) -> SeqDb {
        let cwd = std::env::temp_dir().join(format!("losat_prepass_cwd_{}", std::process::id()));
        std::fs::create_dir_all(&cwd).unwrap();
        let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None).unwrap();
        let path = root().join("fixtures/outfmt16/db").join(name);
        SeqDb::open_in(path.to_str().unwrap(), SeqType::Protein, &search).unwrap()
    }

    fn variant(name: &str) -> PathBuf {
        root().join("seqdb/taxdb").join(name)
    }

    #[test]
    fn names_that_can_be_written_do_not_hold() {
        for name in ["sprot_tax", "sprot_mixtax"] {
            let db = database(name, Some(&variant("good")));
            assert!(db.taxdb_found());
            assert!(!taxdb_may_refuse(&db), "{name}");
        }
        // 10090 has a common name and a kingdom with bytes >= 0x80 in this taxonomy database,
        // which XML2 does not write: only the scientific name is checked.
        let db = database("sprot_mixtax", Some(&variant("high_bytes_latin1")));
        assert!(!taxdb_may_refuse(&db));
    }

    #[test]
    fn a_name_that_cannot_be_written_holds_whether_or_not_the_taxid_is_reported() {
        // sprot_tax: 9606 is the taxid of both sequences. sprot_mixtax: 9606 belongs to the
        // third sequence only, which a search for the first two never reports.
        for name in ["sprot_tax", "sprot_mixtax"] {
            let db = database(name, Some(&variant("high_bytes_utf8")));
            assert!(taxdb_may_refuse(&db), "{name}");
        }
    }

    #[test]
    fn a_record_past_the_last_page_holds() {
        // The record of 10090 (entry 3) ends at the offset of entry 4: 5000 bytes, past the
        // 4096-byte page that ends the mapping of the small `taxdb.btd`.
        let dir = std::env::temp_dir().join(format!("losat_prepass_btd_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::copy(variant("good").join("taxdb.btd"), dir.join("taxdb.btd")).unwrap();
        let mut index = std::fs::read(variant("good").join("taxdb.bti")).unwrap();
        let at = 24 + 8 * 4 + 4;
        index[at..at + 4].copy_from_slice(&5000u32.to_be_bytes());
        std::fs::write(dir.join("taxdb.bti"), &index).unwrap();
        assert!(taxdb_may_refuse(&database("sprot_mixtax", Some(&dir))));
        // 9606 (entry 2) is fine: the database with that taxid alone does not hold.
        assert!(!taxdb_may_refuse(&database("sprot_tax", Some(&dir))));
    }

    #[test]
    fn a_missing_or_damaged_taxdb_that_ncbi_survives_does_not_hold() {
        // No `taxdb.bti` on the search path: the formats write no name.
        let db = database("sprot_mixtax", None);
        assert!(!db.taxdb_found());
        assert!(!taxdb_may_refuse(&db));
        // NCBI turns a wrong magic number into a message and N/A (`expected.tsv`).
        let db = database("sprot_mixtax", Some(&variant("bad_magic")));
        assert!(!taxdb_may_refuse(&db));
    }

    #[test]
    fn the_two_sides_agree_and_the_taxonomy_side_alone_clears_a_good_taxdb() {
        use crate::seqdb::taxdb::LazyTaxDb;
        for (taxonomy, expected) in [
            ("good", Some(false)),
            ("high_bytes_latin1", Some(false)),
            ("high_bytes_utf8", Some(true)),
            ("offset_next_past_end", Some(true)),
            ("bad_magic", Some(false)),
            // Not checked: whatever the pass says, the answer is the database side's.
            ("unsorted_adjacent", None),
        ] {
            for name in ["sprot", "sprot_tax", "sprot_mixtax"] {
                let db = database(name, Some(&variant(taxonomy)));
                let taxdb = LazyTaxDb::new(db.search_path().unwrap().clone());
                let side = taxdb_side(&taxdb);
                if let Some(side) = side {
                    // The taxonomy side is a superset of the database side.
                    assert!(side || !database_side(&db, &taxdb), "{taxonomy} {name}");
                    if !side {
                        assert!(!taxdb_may_refuse(&db), "{taxonomy} {name}");
                    }
                }
                if name == "sprot_mixtax" && expected.is_some() {
                    assert_eq!(side, expected, "{taxonomy}");
                }
                // The database side alone and the dispatch give one answer.
                assert_eq!(
                    taxdb_may_refuse(&db),
                    database_side(&db, &LazyTaxDb::new(db.search_path().unwrap().clone())),
                    "{taxonomy} {name}"
                );
            }
        }
    }

    /// Times the taxonomy side on a real NCBI `taxdb` (`LOSAT_TEST_NCBI_TAXDB=<dir with taxdb.bti
    /// and taxdb.btd>`): `cargo test --release`-like timing of one pass, and whether any record
    /// of it could be refused. Ignored: it needs the 190 MB download of OF5b's fix.
    #[test]
    #[ignore]
    fn real_ncbi_taxdb_walk() {
        use crate::seqdb::taxdb::LazyTaxDb;
        let Some(dir) = std::env::var_os("LOSAT_TEST_NCBI_TAXDB") else {
            return;
        };
        let cwd = std::env::temp_dir();
        let search = SearchPath::new(&cwd, Some(&dir), None).unwrap();
        let taxdb = LazyTaxDb::new(search);
        let entries = taxdb.entry_count().unwrap();
        for _ in 0..3 {
            let start = std::time::Instant::now();
            let answer = taxdb_side(&taxdb);
            eprintln!(
                "real taxdb: {entries} entries, taxonomy side {answer:?} in {:?}",
                start.elapsed()
            );
            assert_eq!(answer, Some(false));
        }
    }

    #[test]
    fn a_database_without_taxids_does_not_hold() {
        let db = database("sprot", Some(&variant("high_bytes_utf8")));
        assert!(db.taxdb_found());
        assert!(!taxdb_may_refuse(&db));
    }
}
