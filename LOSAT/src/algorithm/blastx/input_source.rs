//! BLASTX's reading of FASTA input: the readers of the other programs
//! (`blastinput::fasta_reader`, NCBI's `CBlastFastaInputSource`) read the records, an adapter
//! gives each its BLASTX form (`from_upper` / `from_shared_ref`), and one batch cutter (`BatchCutter`) makes the
//! reader batches of the command line (`QueryStream`) and of the library (`local::reader_batches`).
//! RPSTBLASTN reads its queries through the same stream.
use anyhow::Result;
use std::io::{Read, Seek};

use super::input::{attach_skipped, FastaRecord, SkippedRecord};
use crate::blastinput::fasta_reader::{
    self, FastaInputSource, ParseErrorCode, ReadError, ReaderConfig,
};
use crate::blastinput::query_batch::BatchCount;
use crate::blastinput::seq_range::{record_interval, RecordInterval, SequenceRange};

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1082-1090
/// BLASTX's record from the shared reader's `FastaInputSource::next_sequence_upper`: the residues
/// are already in upper case and the lowercase mask intervals (only when `-lcase_masking` is set)
/// come directly from the reader without a lowercase-and-rescan round trip; the local id, title
/// bytes, and source are moved.
pub(crate) fn from_upper(
    record: fasta_reader::FastaRecord,
    lowercase_masks: Vec<(i32, i32)>,
) -> FastaRecord {
    FastaRecord {
        internal_id: record.local_id,
        title: record.title,
        sequence: record.sequence,
        lowercase_masks,
        source: record.source,
        around: None,
        skipped: None,
    }
}

/// `from_shared` of a record that the caller keeps (the library's records, which the reader of
/// the other programs read or a caller made): its letters are copied, and a record of more than
/// `Int4` letters is refused, as the command line's reader refuses it.
pub(crate) fn from_shared_ref(
    record: &fasta_reader::FastaRecord,
    lowercase: bool,
) -> Result<FastaRecord> {
    if record.sequence.len() > i32::MAX as usize {
        anyhow::bail!("BLASTX sequence exceeds Int4 length");
    }
    Ok(adapt(
        record.local_id.clone(),
        record.title.clone(),
        record.sequence.clone(),
        record.source.clone(),
        lowercase,
    ))
}

fn adapt(
    internal_id: String,
    title: Vec<u8>,
    mut sequence: Vec<u8>,
    source: Option<Box<fasta_reader::RecordSource>>,
    lowercase: bool,
) -> FastaRecord {
    let mut lowercase_masks: Vec<(i32, i32)> = Vec::new();
    // A record without lower case (nearly every one) is a read of its letters and no more.
    if sequence.iter().any(u8::is_ascii_lowercase) {
        let mut next = 0;
        while let Some(start) = sequence[next..]
            .iter()
            .position(u8::is_ascii_lowercase)
            .map(|at| next + at)
        {
            let end = sequence[start..]
                .iter()
                .position(|letter| !letter.is_ascii_lowercase())
                .map_or(sequence.len(), |at| start + at);
            sequence[start..end].make_ascii_uppercase();
            if lowercase {
                // The length is at most `i32::MAX` (the readers refuse more).
                lowercase_masks.push((start as i32, end as i32 - 1));
            }
            next = end;
        }
    }
    FastaRecord {
        internal_id,
        title,
        sequence,
        lowercase_masks,
        source,
        around: None,
        skipped: None,
    }
}

/// How a reader batch is cut: the range that cuts each record (`-query_loc`), whether a cut
/// record keeps its whole letters (the archive, `FastaRecord::cut_in_place`), the letters on each
/// side of its interval that it keeps (`FastaRecord::keep_around`), whether the batch keeps the
/// records that the range skips (`FastaRecord::skipped`), the letters of a batch, and under
/// `-parse_deflines` NCBI's count of the records (`BatchCount`).
pub(crate) struct CutSettings<'r> {
    pub(crate) range: Option<&'r SequenceRange>,
    pub(crate) keep_whole: bool,
    pub(crate) around: usize,
    pub(crate) keep_skipped: bool,
    pub(crate) batch_size: usize,
    pub(crate) count: Option<BatchCount>,
    /// Whether an input of no record (blank lines alone) gives one empty batch (NCBI's query
    /// factory refuses it, `Empty CBlastQueryVector`; every caller but ABI v1's).
    pub(crate) blank_batch: bool,
}

type Deliver<'d> = &'d mut dyn FnMut(&[FastaRecord]) -> Result<()>;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-171
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-293
/// The reader batches of BLASTX's records (`add` each record in the order it was read, `finish`
/// at the end), given to `deliver` in order: the batch once it holds `batch_size` letters or more,
/// the rest at the end. With a range (`-query_loc`) the batch counts the length of the whole
/// record and cuts the record to its interval; a record whose interval starts more than one
/// letter past its end adds nothing and is not kept (it keeps its number: the reader counted it),
/// so the last batch can hold none (NCBI's query factory refuses it, `Empty CBlastQueryVector`).
/// Under `-parse_deflines` a record counts as NCBI's reader counts it (`BatchCount`: a repeated
/// id by its first record's letters). Each delivered batch is dropped: one batch is held at a time,
/// as NCBI's loop holds one.
pub(crate) struct BatchCutter<'r> {
    settings: CutSettings<'r>,
    batch: Vec<FastaRecord>,
    /// The batch's records that the range skipped, for its first record (`FastaRecord::skipped`):
    /// NCBI's reader added them to its scope.
    skipped: Vec<SkippedRecord>,
    /// The letters the batch has read (`size_read`).
    letters: usize,
    /// The records the batch has read, kept or not.
    read: usize,
    /// Whether a batch has been delivered.
    delivered: bool,
}

impl<'r> BatchCutter<'r> {
    pub(crate) fn new(settings: CutSettings<'r>) -> Self {
        Self {
            settings,
            batch: Vec::new(),
            skipped: Vec::new(),
            letters: 0,
            read: 0,
            delivered: false,
        }
    }

    /// The next record read: cut to its interval and added to the batch, which is delivered when
    /// it holds `batch_size` letters.
    pub(crate) fn add(&mut self, mut record: FastaRecord, deliver: Deliver<'_>) -> Result<()> {
        let settings = &mut self.settings;
        let length = record.sequence.len();
        self.read += 1;
        let interval = settings.range.map(|range| record_interval(range, length));
        let skipped = matches!(interval, Some(RecordInterval::PastEnd));
        // A skipped record counts 0 and joins the scope (`BatchCount::count`).
        self.letters += match settings.count.as_mut() {
            Some(count) => count.count(record.believed(), length, skipped),
            None if skipped => 0,
            None => length,
        };
        match interval {
            None => self.batch.push(record),
            Some(RecordInterval::Letters { from, to }) => {
                // The letter before the interval is the only one outside it that a split query's
                // chunks read (`split::x_chunk_masks`); a record keeps its whole letters for it.
                let keep_whole = settings.keep_whole || record.lower_case_before(from);
                record.keep_around(from, to + 1, settings.around);
                record.cut_in_place(from, to + 1, keep_whole);
                self.batch.push(record);
            }
            Some(RecordInterval::Empty) => {
                record.cut_in_place(length, length, settings.keep_whole);
                self.batch.push(record);
            }
            Some(RecordInterval::PastEnd) => {
                if settings.keep_skipped {
                    self.skipped
                        .push(SkippedRecord::of(&record, self.batch.len()));
                }
            }
        }
        if self.letters >= self.settings.batch_size {
            self.deliver_batch(deliver)?;
            self.letters = 0;
            self.read = 0;
            if let Some(count) = self.settings.count.as_mut() {
                count.end_batch();
            }
        }
        Ok(())
    }

    /// The end of the records: the batch being read is delivered, even when it holds none because
    /// the range skipped its records (NCBI's query factory then refuses it). Input of no record
    /// (blank lines alone) gives one empty batch; the caller does not call this for an input that
    /// has nothing at all.
    pub(crate) fn finish(&mut self, deliver: Deliver<'_>) -> Result<()> {
        if self.read > 0 || (self.settings.blank_batch && !self.delivered) {
            self.deliver_batch(deliver)?;
        }
        Ok(())
    }

    fn deliver_batch(&mut self, deliver: Deliver<'_>) -> Result<()> {
        attach_skipped(&mut self.batch, &mut self.skipped);
        self.delivered = true;
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-261
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:931-938
        // NCBI's batch is local to one turn of the loop: the batch read for the next turn
        // replaces it.
        let result = deliver(&self.batch);
        self.batch.clear();
        result
    }
}

/// How the search reads its query file: `-parse_deflines`, the search's ids and `-strand`
/// (`fasta_reader::ReaderConfig`). The read of the file before the search
/// (`blast_xml::prescan_queries`) reads it the same way, so that it sees the search's records and
/// titles and stops at the reader error that stops the search (under `-parse_deflines`, a first
/// line that is not a defline).
pub(crate) struct QueryStreamSettings {
    pub(crate) parse_deflines: bool,
    pub(crate) id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    pub(crate) strand: crate::report::archive::bioseq::SeqLocStrand,
    /// NCBI's data loaders (`ReaderConfig::data_loaders`).
    pub(crate) data_loaders: bool,
    /// The program of the reader's messages (`args::ProgramName`; BLASTX unless the run is
    /// RPSTBLASTN's).
    pub(crate) program: super::args::ProgramName,
    /// NCBI's count of the records into the batches (`QueryStream::batch_count`), from the first
    /// record on.
    pub(crate) batch_count: Option<BatchCount>,
}

impl QueryStreamSettings {
    /// The search's stream of the query file; `keep_defline`: the BLAST archive (outfmt 11) or
    /// the search strategy writes each record's defline line, and a record cut to a range keeps
    /// its whole letters.
    pub(crate) fn stream(
        &self,
        file: impl Into<crate::blastinput::input_files::InputFile>,
        keep_defline: bool,
    ) -> QueryStream<crate::blastinput::input_files::InputFile> {
        let mut config = ReaderConfig::query(self.program.refusal, false, self.data_loaders)
            .with_defline(keep_defline)
            .with_parse_deflines(self.parse_deflines);
        config.strand = self.strand;
        let mut source = FastaInputSource::from_file(file, config);
        source.set_id_scope(self.id_scope.clone());
        QueryStream {
            source,
            keep_whole: keep_defline,
            around: 0,
            // The records that `-query_loc` skips join the scopes of the queries
            // (`input::scope_records`).
            keep_skipped: self.parse_deflines,
            batch_count: self.batch_count.clone(),
        }
    }

    /// A stream of the query file for a read before the search: its ids go to a copy of the
    /// search's (`IdScope::copy`), so that it stops at the record the search refuses (a query
    /// with the id of a subject) and the search reads its records into its ids as before.
    pub(crate) fn stream_before_search(
        &self,
        file: impl Into<crate::blastinput::input_files::InputFile>,
    ) -> QueryStream<crate::blastinput::input_files::InputFile> {
        let mut stream = self.stream(file, false);
        stream.source.set_id_scope(
            self.id_scope
                .as_deref()
                .map(|scope| std::sync::Arc::new(scope.copy())),
        );
        stream
    }
}

/// The query records of one input, read one at a time and cut into reader batches
/// (`read_batches`).
pub(crate) struct QueryStream<R: Read> {
    source: FastaInputSource<R>,
    /// Whether a cut record keeps its whole letters (`CutSettings::keep_whole`).
    keep_whole: bool,
    /// The letters on each side of its `-query_loc` interval that each record cut to it keeps
    /// (`FastaRecord::keep_around`, `empty_rows::query_letters_around`); 0 for none.
    pub(crate) around: usize,
    /// Whether the reader keeps the records that the `-query_loc` range skips
    /// (`FastaRecord::skipped`): NCBI's scope of a batch's prefetch holds them, which only the
    /// stop of a prefetch asks (`common::database::PrefetchQueries`), so the caller sets it where
    /// a prefetch of the run can stop (`BlastxDatabase::prefetch_can_stop`; review 5 of line DBID
    /// loadord, L2).
    pub(crate) keep_skipped: bool,
    /// Under `-parse_deflines`, NCBI's count of the records into the reader batches
    /// (`BatchCount`); `None` counts each record by its own letters.
    batch_count: Option<BatchCount>,
}

impl<R: Read + Seek> QueryStream<R> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
    /// Whether the input holds white space only (`FastaInputSource::stream_is_empty`), before the
    /// first record is read.
    pub(crate) fn stream_is_empty(&mut self) -> bool {
        self.source.stream_is_empty()
    }
}

impl<R: Read> QueryStream<R> {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    pub(crate) fn at_eof(&mut self) -> bool {
        self.source.end()
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873 (cited at the read
    // before the search, `blast_xml::prescan_queries`: NCBI reads ahead on the argument's one stream
    // and puts it back at `orig_p`)
    /// The input this stream reads, given back: the bytes it buffered are dropped, so the caller
    /// puts the input's position where its next reading starts (`blast_xml::prescan_queries`, which
    /// reads the query file before the search and hands the same open file to the search).
    pub(crate) fn into_input(self) -> R {
        self.source.into_input()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:134-171
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-262
    /// Reads the records of the input to its end, cut into the reader batches of `batch_size`
    /// letters (`BatchCutter`), given to `deliver` in order, each record's messages going to
    /// `warn` as the record is read. Each batch is dropped once delivered, so that the read holds
    /// one batch, not the file (as NCBI's loop holds one; the read of the queries before the
    /// search, `blast_xml::prescan_queries`, has no NCBI counterpart). A reader error ends the
    /// read in the batch being read, which is never delivered, after the messages of its records
    /// and of the record the reader failed on; blank lines alone, an input without a record, are
    /// one batch without a query (NCBI's `Empty CBlastQueryVector`); an input without any line
    /// gives no batch.
    pub(crate) fn read_batches(
        &mut self,
        lowercase: bool,
        batch_size: usize,
        range: Option<&SequenceRange>,
        deliver: Deliver<'_>,
        warn: &mut dyn FnMut(&[u8]) -> Result<()>,
    ) -> Result<()> {
        let settings = CutSettings {
            range,
            keep_whole: self.keep_whole,
            around: self.around,
            keep_skipped: self.keep_skipped,
            batch_size,
            count: self.batch_count.take(),
            blank_batch: true,
        };
        read_source_batches(
            &mut self.source,
            lowercase,
            settings,
            deliver,
            warn,
            ReadError::into_app_error,
        )
    }
}

/// The records of `source`, read one at a time and cut into reader batches (`BatchCutter`), given
/// to `deliver`; an input without a line gives no batch. A reader error is turned into the run's
/// error by `error` (the command line's: `ReadError::into_app_error`; ABI v1's: the reader's
/// message as it is).
pub(crate) fn read_source_batches<R: Read>(
    source: &mut FastaInputSource<R>,
    lowercase: bool,
    settings: CutSettings<'_>,
    deliver: Deliver<'_>,
    warn: &mut dyn FnMut(&[u8]) -> Result<()>,
    error: fn(ReadError) -> anyhow::Error,
) -> Result<()> {
    if source.end() {
        return Ok(());
    }
    let mut cutter = BatchCutter::new(settings);
    while !source.end() {
        match next_record(source, lowercase, warn) {
            Ok(Some(record)) => cutter.add(record, deliver)?,
            Ok(None) => break,
            Err(Failure::Write(written)) => return Err(written),
            Err(Failure::Reader(failure)) => return Err(error(failure)),
        }
    }
    cutter.finish(deliver)
}

/// Why `next_record` read no record.
enum Failure {
    /// The reader's error.
    Reader(ReadError),
    /// The error of the writer of the messages.
    Write(anyhow::Error),
}

/// The next record of `source` with its messages written to `warn`; `None` at the end the reader
/// signals with `eEOF` (only blank and comment lines were left). An `Err` is the reader's error,
/// or the error of `warn` where the messages could not be written.
fn next_record<R: Read>(
    source: &mut FastaInputSource<R>,
    lowercase: bool,
    warn: &mut dyn FnMut(&[u8]) -> Result<()>,
) -> std::result::Result<Option<FastaRecord>, Failure> {
    let mut written = None;
    let read = source.next_sequence_upper(lowercase, &mut |message| {
        warn(message).map_err(|error| {
            written = Some(error);
            std::io::Error::other("a reader message was not written")
        })
    });
    match read {
        Ok((record, masks)) => Ok(Some(from_upper(record, masks))),
        Err(_) if written.is_some() => Err(Failure::Write(written.take().expect("checked"))),
        Err(ReadError::Parse {
            code: ParseErrorCode::Eof,
            ..
        }) => Ok(None),
        Err(error) => Err(Failure::Reader(error)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:302-314
/// The records of LOSAT Web's ABI v1 input `data` (a subject), as the readers of the other programs
/// read bytes in memory (`FastaInputSource::from_bytes`), without the data loaders. A reader error
/// is the reader's message as it is (ABI v1 is frozen).
#[cfg(any(losat_web_abi_v1, test))]
pub(crate) fn read_web_subjects(
    data: &[u8],
    lowercase: bool,
    warn: &mut dyn FnMut(&[u8]) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    let mut source =
        FastaInputSource::from_bytes(data, ReaderConfig::subject("BLASTX", true, false));
    let mut records = Vec::new();
    while !source.end() {
        match next_record(&mut source, lowercase, warn) {
            Ok(Some(record)) => records.push(record),
            Ok(None) => break,
            Err(Failure::Write(written)) => return Err(written),
            Err(Failure::Reader(failure)) => return Err(failure.into()),
        }
    }
    Ok(records)
}

/// The query batches of LOSAT Web's ABI v1 input `data` (`read_source_batches`): no data loaders, the
/// reader's message as it is for an error, and blank lines alone give no batch (ABI v1 is frozen).
#[cfg(any(losat_web_abi_v1, test))]
pub(crate) fn read_web_queries(
    data: &[u8],
    lowercase: bool,
    batch_size: usize,
    deliver: Deliver<'_>,
    warn: &mut dyn FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    let mut source =
        FastaInputSource::from_bytes(data, ReaderConfig::query("BLASTX", false, false));
    let settings = CutSettings {
        range: None,
        keep_whole: false,
        around: 0,
        keep_skipped: false,
        batch_size,
        count: None,
        blank_batch: false,
    };
    read_source_batches(
        &mut source,
        lowercase,
        settings,
        deliver,
        warn,
        anyhow::Error::from,
    )
}

#[cfg(test)]
pub(crate) use tests::{parse_fasta, read_fasta};

#[cfg(test)]
#[path = "input_source_tests.rs"]
mod tests;
