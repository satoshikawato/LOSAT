//! NCBI's sequence streams of index creation (`CSequenceIStream`): the sequences of a BLAST
//! database (`CSequenceIStreamBlastDB`, with the masks of one filtering algorithm) and of a
//! FASTA file (`CSequenceIStreamFasta`, with the lower case letters as masks). Each `next`
//! gives one sequence or the end of the input; `putback` makes the next `next` give the last
//! sequence again.

use std::io::Read;
use std::rc::Rc;

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::seqdb::mask::MaskRanges;
use crate::seqdb::MaskBlock;
use crate::seqdb::SeqType;
use crate::seqid::model::ParsedSeqId;
use crate::seqmasks_io::bdb_reader::{BdbError, MaskBdbReader};
use crate::sequence::title::get_title_without_source;
use crate::writedb::fasta_source::reader::{
    past_seq_pos_refusal, FastaReader, ReadError, ReaderFlags, ReaderWarning, Record, SeqInst,
};

/// `CSequenceIStream::TSeqData` of a sequence (a null `seq_entry_` is `None` of `next`).
#[derive(Debug, Default)]
pub(crate) struct SeqData {
    /// The IUPACna letters of the sequence (its `CSeqVector` in `eCoding_Iupac`), in the form
    /// the index reads them.
    pub(crate) bases: Bases,
    /// `mask_locs_` that are packed intervals (the only ones `CMaskHelper::Add` keeps), in their
    /// order: each one's intervals `(from, to)`, `to` inclusive.
    pub(crate) masks: Vec<Vec<(u32, u32)>>,
    /// The first word of `sequence::GetTitle` of the Bioseq handle (`extractSeqVector`), when the
    /// stream was asked for titles (`-idmap`; nothing else reads it).
    pub(crate) id: Option<Vec<u8>>,
    /// The messages that NCBI's `GetBioseq` posts while it makes the sequence (taxonomy
    /// lookups), to be posted before it is used.
    pub(crate) texts: Vec<String>,
}

/// The bit of `CODES` that marks a letter other than A, C, G or T.
const INVALID: u8 = 4;

/// One less than `base_value` of each IUPACna letter `r` as the index packs it: the NCBI2na value
/// of `A`, `C`, `G`, `T`, and 0 with `INVALID` for any other letter (`base_value` 0).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:103-112
const CODES: [u8; 256] = {
    let mut codes = [INVALID; 256];
    codes[b'A' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'G' as usize] = 2;
    codes[b'T' as usize] = 3;
    codes
};

/// The letters of a byte of `Bases::packed`.
const CR: usize = super::CR as usize;

/// The letters of a sequence as the index reads them: whether each is A, C, G or T
/// (NCBI's `base_value` above 0), and the sequence packed as `seq_store_` holds it. The letters
/// themselves (one byte each) are not kept.
#[derive(Debug, Default)]
pub(crate) struct Bases {
    /// The number of letters.
    len: usize,
    /// Four letters to a byte, the first in the two highest bits, each `CODES[r] & 3` (0 for a
    /// letter other than A, C, G, T), the last byte's unused bits 0.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:814-876
    packed: Vec<u8>,
    /// The runs `[from, to)` of letters other than A, C, G and T, in increasing order.
    invalid: Vec<(u32, u32)>,
    /// The letters pushed into the last byte not yet in `packed`, and their number.
    accum: u8,
    pending: u8,
}

impl Bases {
    /// The bases of `letters`.
    pub(crate) fn of_letters(letters: &[u8]) -> Self {
        let mut bases = Bases::default();
        bases.push(letters);
        bases.finish();
        bases
    }

    /// Adds the next letters of the sequence: the letters that complete a partly filled byte one
    /// at a time, then each whole group of four as one byte, then the rest one at a time.
    fn push(&mut self, letters: &[u8]) {
        let mut rest = letters;
        while self.pending != 0 {
            let Some((&r, after)) = rest.split_first() else {
                return;
            };
            self.push_letter(r);
            rest = after;
        }
        let (groups, tail) = rest.as_chunks::<CR>();
        let start = self.len;
        let invalid = &mut self.invalid;
        // One byte more for the letters that `finish` packs.
        self.packed.reserve(groups.len() + 1);
        self.packed
            .extend(groups.iter().enumerate().map(|(group, letters)| {
                let codes = letters.map(|r| CODES[usize::from(r)]);
                if (codes[0] | codes[1] | codes[2] | codes[3]) & INVALID != 0 {
                    for (k, &code) in codes.iter().enumerate() {
                        if code & INVALID != 0 {
                            add_invalid(invalid, (start + group * CR + k) as u32);
                        }
                    }
                }
                codes.iter().fold(0, |byte, &code| (byte << 2) | (code & 3))
            }));
        self.len += groups.len() * CR;
        for &r in tail {
            self.push_letter(r);
        }
    }

    /// Adds one letter.
    fn push_letter(&mut self, r: u8) {
        let code = CODES[usize::from(r)];
        if code & INVALID != 0 {
            add_invalid(&mut self.invalid, self.len as u32);
        }
        self.accum = (self.accum << 2) | (code & 3);
        self.pending += 1;
        if usize::from(self.pending) == CR {
            self.packed.push(self.accum);
            (self.accum, self.pending) = (0, 0);
        }
        self.len += 1;
    }

    /// Ends the sequence: its last letters fill the high bits of the last byte.
    fn finish(&mut self) {
        if self.pending > 0 {
            self.packed
                .push(self.accum << ((CR as u8 - self.pending) * 2));
            (self.accum, self.pending) = (0, 0);
        }
        self.packed.shrink_to_fit();
        self.invalid.shrink_to_fit();
    }

    /// The number of letters (`c_seq_.size()` before NCBI's `TSeqPos` cast).
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// The packed letters (`seq_store_`'s bytes of the sequence).
    pub(crate) fn packed(&self) -> &[u8] {
        &self.packed
    }

    /// The runs of letters other than A, C, G and T (NCBI's `base_value` 0).
    pub(crate) fn invalid(&self) -> &[(u32, u32)] {
        &self.invalid
    }
}

/// Adds the letter at `pos` to the runs of letters other than A, C, G and T.
fn add_invalid(invalid: &mut Vec<(u32, u32)>, pos: u32) {
    match invalid.last_mut() {
        Some(run) if run.1 == pos => run.1 = pos + 1,
        _ => invalid.push((pos, pos + 1)),
    }
}

/// The letters a database sequence is read by at a time, so that it is never held whole as
/// letters.
const READ_BLOCK: usize = 1 << 20;

/// Why a stream stopped.
#[derive(Debug)]
pub(crate) enum StreamError {
    /// An exception of NCBI's, which the application reports (`Application's execution failed`).
    Ncbi(NcbiException),
    /// LOSAT's refusal, the whole line after `Error: `.
    Refused(String),
}

/// `CSequenceIStream`.
pub(crate) trait SequenceStream {
    /// `next()`: the next sequence, `None` at the end of the input.
    fn next(&mut self) -> Result<Option<Rc<SeqData>>, StreamError>;
    /// `putback()`.
    fn putback(&mut self);
    /// The lines the reader posted (`LOG_POST`) since the last call, to be posted now.
    fn take_warnings(&mut self) -> Vec<Vec<u8>> {
        Vec::new()
    }
}

/// `CSequenceIStream::CSequenceIStream_Exception::eParam` (`database parameter error`).
const PARAM_ERROR: &str = "CSequenceIStream_Exception::database parameter error";

/// The first word of a title: what comes before its first space or tab.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:658-661
fn first_word(title: &[u8]) -> Vec<u8> {
    let end = title
        .iter()
        .position(|&byte| byte == b' ' || byte == b'\t')
        .unwrap_or(title.len());
    title[..end].to_vec()
}

/// `CSequenceIStreamBlastDB`: the sequences of a BLAST database, from its OID 0
/// (`CSeqDB::GetBioseq`), each with the masks of one filtering algorithm when one is used.
pub(crate) struct BlastDbStream {
    reader: MaskBdbReader,
    /// `filter_algo_id_` when `use_filter_`.
    filter: Option<i32>,
    /// The buffers of the masks' reads, kept from one sequence to the next
    /// (`SeqDb::mask_ranges_in`).
    mask_block: MaskBlock,
    /// NCBI's warning 803.10 of the read of a header that took a UTF8String for a VisibleString
    /// (without its line feed), once in the run, to be posted before the sequence's messages.
    warnings: Vec<Vec<u8>>,
}

/// The exception or refusal of a database error of the reader.
fn bdb_error(error: BdbError) -> StreamError {
    match error {
        BdbError::Ncbi(exception) => StreamError::Ncbi(exception),
        BdbError::Refused(reason) => StreamError::Refused(reason),
        BdbError::Unsupported(what) => StreamError::Refused(super::unsupported(&what)),
    }
}

/// How a filtering algorithm is named: by its id (the `CSequenceIStreamBlastDB` constructor with
/// an `int`) or by the text of `-db_mask` (the `string` constructor).
pub(crate) enum Filter<'a> {
    None,
    /// The bytes of `-db_mask`.
    Name(&'a [u8]),
}

impl BlastDbStream {
    /// `CSequenceIStreamBlastDB` with `dbname`, `use_filter` and `filter`: opens the database
    /// (`CSeqDB` with `dbname` and `eNucleotide`), then checks the filtering algorithm
    /// (`CheckBlastDBMaskAlgorithmId`), or finds the id of its name (`GetMaskAlgorithmId`): a name
    /// that is a number (`NStr::StringToNumeric`, base 10) is an id.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_bdb.cpp:48-57
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_bdb.cpp:62-96
    pub(crate) fn open(
        dbname: &[u8],
        filter: Filter<'_>,
        registry: Option<&[u8]>,
        titles: bool,
    ) -> Result<Self, StreamError> {
        let mut reader =
            MaskBdbReader::open(dbname, SeqType::Nucleotide, registry, false, super::PROGRAM)
                .map_err(bdb_error)?;
        if titles {
            reader.ask_get_title();
        }
        let mut stream = BlastDbStream {
            reader,
            filter: None,
            mask_block: MaskBlock::default(),
            warnings: Vec::new(),
        };
        stream.filter = match filter {
            Filter::None => None,
            Filter::Name(name) => match crate::ncbi_str::string_to_int(name, 0).ok() {
                Some(id) => Some(stream.checked_id(id)?),
                None => Some(stream.named_id(name)?),
            },
        };
        Ok(stream)
    }

    /// `FindVolumePaths` (with `paths` and true) of the database opened as `dbname`.
    pub(crate) fn volume_paths(&self, dbname: &[u8]) -> Result<Vec<Vec<u8>>, StreamError> {
        self.reader.volume_paths(dbname).map_err(bdb_error)
    }

    /// `GetNumOIDs()` of the database.
    pub(crate) fn num_oids(&self) -> u32 {
        self.reader.db().num_oids()
    }

    /// `GetAvailableMaskAlgorithmDescriptions()` of the database: the descriptions of its
    /// filtering algorithms, empty when it has none (also `ShowSupportedFilters` of `dbname`).
    // NCBI reference (598d8ae6): c++/include/algo/blast/dbindex/sequence_istream_bdb.hpp:53-57
    pub(crate) fn descriptions(&self) -> Result<String, StreamError> {
        self.reader
            .db()
            .mask_algorithms()
            .and_then(|algorithms| algorithms.descriptions())
            .map_err(|error| bdb_error(self.reader.database_error(error)))
    }

    /// `CheckBlastDBMaskAlgorithmId` with `db` and `id`: `id` when the database has that algorithm.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_bdb.cpp:48-57
    fn checked_id(&self, id: i32) -> Result<i32, StreamError> {
        let ids = self
            .reader
            .db()
            .mask_algorithms()
            .map(|algorithms| algorithms.ids())
            .map_err(|error| bdb_error(self.reader.database_error(error)))?;
        if ids.contains(&id) {
            return Ok(id);
        }
        Err(StreamError::Ncbi(NcbiException::new(
            ThrowSite {
                file: "src/algo/blast/dbindex/sequence_istream_bdb.cpp",
                line: 53,
                kind: PARAM_ERROR,
                function: None,
            },
            format!("unrecognized filter algorithm id{}", self.descriptions()?),
        )))
    }

    /// `GetMaskAlgorithmId` of the database (`seqdb_`) with `name`; its `CSeqDBException` is thrown
    /// again as a `CSequenceIStream_Exception`. Its text is the reader's (`MaskAlgorithms::id_of`):
    /// the remapper's for mask columns, the names table of `MASKLIST` for GI-indexed masks.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_bdb.cpp:78-94
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2261-2266
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2292-2296
    fn named_id(&self, name: &[u8]) -> Result<i32, StreamError> {
        let found = self
            .reader
            .db()
            .mask_algorithms()
            .map(|algorithms| algorithms.id_of(name))
            .map_err(|error| bdb_error(self.reader.database_error(error)))?;
        let message = match found {
            Ok(id) => return Ok(id),
            Err(message) => message,
        };
        let site = crate::seqdb::throw_site(&message)
            .unwrap_or(ThrowSite::generated("CSeqDBException::eInvalid"));
        Err(StreamError::Ncbi(
            NcbiException::new(site, message).rethrow(
                ThrowSite {
                    file: "src/algo/blast/dbindex/sequence_istream_bdb.cpp",
                    line: 90,
                    kind: PARAM_ERROR,
                    function: None,
                },
                format!("unrecognised filter algorithm name{}", self.descriptions()?),
            ),
        ))
    }
}

impl BlastDbStream {
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:296-310
    /// NCBI's warning 803.10 of the read of the header of `oid` (`GetBioseq`), when it met a
    /// UTF8String for a VisibleString before any failure and the run has not posted it
    /// (`SeqDb::take_utf8_post_at`): without a prefix, as makembindex's diagnostics print it
    /// (NCBI 2.17.0 probe of line RG item T5d, merge 18).
    #[cold]
    fn post_utf8_string(&mut self, oid: u32) {
        let db = self.reader.db();
        let Some(path) = db.location_of(oid).and_then(|at| db.take_utf8_post_at(at)) else {
            return;
        };
        let line = crate::seqid::ber::utf8_string_warning(None, &path);
        self.warnings
            .push(line.strip_suffix(b"\n").unwrap_or(&line).to_vec());
    }
}

impl SequenceStream for BlastDbStream {
    /// `CSequenceIStreamBlastDB::next()`: the sequence of the current OID and, with a filter,
    /// its masks as one packed location (`GetMaskData`: `[first, second)` ranges become
    /// `first..=second - 1`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_bdb.cpp:99-135
    fn next(&mut self) -> Result<Option<Rc<SeqData>>, StreamError> {
        let oid = self.reader.next_oid();
        let read = self.reader.next_record_by_ranges();
        if self.reader.take_utf8_header() {
            self.post_utf8_string(oid);
        }
        let Some((record, ambiguity)) = read.map_err(bdb_error)? else {
            return Ok(None);
        };
        let length = ambiguity.length();
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1358-1386
        // `GetBioseq` gives a null Bioseq for a sequence without letters, and for one whose every
        // definition line the alias files' `OID_MASK_TYPE` erased (the reader's record without
        // letters, line RG D1b-2: probe mk_old_ntall), which `next` dereferences: a
        // `CCoreException` (`eNullPtr`) whose report prints a stack trace with the run's memory
        // addresses (PD-LOSAT-NCBI-DEFECTS rule 3).
        if length == 0 {
            return Err(StreamError::Refused(super::unsupported(&format!(
                "sequence {oid} of the database, which has no letters or no definition line left by \
                 the alias files' OID_MASK_TYPE (NCBI's makembindex reads the null sequence its \
                 database reader gives for it and stops with a report of the run's memory \
                 addresses),"
            ))));
        }
        let mut bases = Bases::default();
        let mut from = 0;
        while from < length {
            let to = length.min(from + READ_BLOCK);
            let letters = self
                .reader
                .db()
                .nucleotide_iupac_range_with(oid, from, to, &ambiguity)
                .map_err(|error| bdb_error(self.reader.database_error(error)))?;
            bases.push(&letters);
            from = to;
        }
        bases.finish();
        let mut masks = Vec::new();
        if let Some(algo_id) = self.filter {
            let db = self.reader.db();
            let ranges = db
                .locate_oid(oid)
                .and_then(|at| db.mask_ranges_in(at, algo_id, &mut self.mask_block))
                .map_err(|error| bdb_error(self.reader.database_error(error)))?;
            match ranges {
                MaskRanges::Ranges(ranges) => masks.push(
                    ranges
                        .into_iter()
                        .map(|(start, end)| (start, end.wrapping_sub(1)))
                        .collect(),
                ),
                MaskRanges::AlgorithmNotInVolume { .. } => {
                    return Err(StreamError::Refused(super::unsupported(&format!(
                        "sequence {oid} of the database, whose volume has no mask data of \
                         algorithm {algo_id} where NCBI reads another volume's,"
                    ))))
                }
            }
        }
        let id = record.get_title.as_deref().map(first_word);
        Ok(Some(Rc::new(SeqData {
            bases,
            masks,
            id,
            texts: record.texts,
        })))
    }

    /// `putback()`: the last OID is read again (none before OID 0).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_bdb.cpp:138-139
    fn putback(&mut self) {
        self.reader.putback();
    }

    fn take_warnings(&mut self) -> Vec<Vec<u8>> {
        std::mem::take(&mut self.warnings)
    }
}

/// `CSequenceIStreamFasta`: the records of a FASTA file or of the standard input, read by
/// `CFastaReader` with makembindex's flags (`ReaderFlags::dbindex_fasta`). An exception of the
/// reader ends the input (`next` catches every exception).
pub(crate) struct FastaStream<'a> {
    reader: FastaReader<Box<dyn Read + 'a>>,
    /// `cache_`: the last sequence `next` gave, and `use_cache_`.
    cache: Option<Rc<SeqData>>,
    use_cache: bool,
    /// `-idmap`: the records' titles are read.
    titles: bool,
    /// The reader's warnings, to be posted (`LOG_POST`, one line each).
    warnings: Vec<Vec<u8>>,
}

impl<'a> FastaStream<'a> {
    /// `CSequenceIStreamFasta` over the opened input (`input`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_fasta.cpp:55-87
    pub(crate) fn new(input: Box<dyn Read + 'a>, titles: bool) -> Self {
        FastaStream {
            reader: FastaReader::new(input, ReaderFlags::dbindex_fasta()),
            cache: None,
            use_cache: false,
            titles,
            warnings: Vec::new(),
        }
    }

    /// One record of the reader with its masks, or `None` when the reader ended or threw.
    fn read(&mut self) -> Result<Option<SeqData>, StreamError> {
        if self.reader.at_eof() {
            let _ = self.reader.take_read_error();
            return Ok(None);
        }
        let warnings = &mut self.warnings;
        let result = self.reader.read_one_seq(&mut |warning| {
            if let ReaderWarning::Line(line) = warning {
                warnings.push(line.to_vec());
            }
        });
        let masks = self.reader.take_masks();
        if self.reader.take_read_error().is_some() {
            return Ok(None);
        }
        let record: Record = match result {
            Ok(record) => record,
            // NCBI catches every exception of the reader and ends the input.
            Err(ReadError::Parse(_) | ReadError::Ncbi(_) | ReadError::OutOfMemory) => {
                return Ok(None)
            }
        };
        if let SeqInst::PastSeqPos { length, .. } = record.inst {
            let mut id = Vec::new();
            if let Some(first) = record.ids.first() {
                crate::seqid::label::fasta(first, &mut id);
            }
            return Err(StreamError::Refused(super::unsupported(
                &past_seq_pos_refusal(&String::from_utf8_lossy(&id), length),
            )));
        }
        let id = if self.titles {
            let ids = [ParsedSeqId::Local(crate::seqid::model::ObjectId::Id(
                record.local_id as i32,
            ))];
            Some(first_word(&get_title_without_source(
                &ids,
                record.title_descr.as_deref().unwrap_or_default(),
            )))
        } else {
            None
        };
        Ok(Some(SeqData {
            bases: Bases::of_letters(&crate::seqmasks_io::fasta_reader::iupac_letters(
                record.inst,
                &mut self.reader,
            )),
            masks: if masks.is_empty() {
                Vec::new()
            } else {
                vec![masks]
            },
            id,
            texts: Vec::new(),
        }))
    }
}

impl SequenceStream for FastaStream<'_> {
    /// `CSequenceIStreamFasta::next()`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_fasta.cpp:98-118
    fn next(&mut self) -> Result<Option<Rc<SeqData>>, StreamError> {
        if self.use_cache {
            self.use_cache = false;
            return Ok(self.cache.clone());
        }
        let result = self.read()?.map(Rc::new);
        self.cache = result.clone();
        Ok(result)
    }

    /// `putback()`: the next `next` gives the last sequence again.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_fasta.cpp:121-122
    fn putback(&mut self) {
        self.use_cache = true;
    }

    fn take_warnings(&mut self) -> Vec<Vec<u8>> {
        std::mem::take(&mut self.warnings)
    }
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
