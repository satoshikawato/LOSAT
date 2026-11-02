//! EXPERIMENT (LOSAT_X_BNQPACK; in the strict set since 2026-10-25: on unless it is `0`): the query
//! records read before the search held in a store that the batch loop takes them out of, their
//! letters packed from the reading on, made again when the loop takes their batch and dropped once
//! the batch has encoded them, where nothing reads them after the search, instead of holding every
//! record's letters until the end of the search. Until 2026-11-01 the store had a switch of its own
//! (LOSAT_X_BNQDROP) and packed only files of two of the task's largest batches or more; both were
//! folded into this switch then (line SF S3, family bn-query-read: below two batches packing read
//! the same time and memory, see `docs/agents/diagnostics.md`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//! 	    BLAST_PROF_START( APP.LOOP.PRE );
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//! ```
//! NCBI holds a batch's Bioseqs in its scope until the batch is formatted, beside the search's
//! own encoded copy (`SetupQueries_OMF`). No NCBI counterpart for the drop: it is memory only.
//!
//! What it changes. Where the query records are read before the search (`search_cli` without
//! LOSAT_X_QSTREAM acting), `search_cli` held every record, with its IUPAC letters, until the
//! search ended, and `search_query_batch` encodes each batch's letters into its BLASTNA buffer (both
//! strands and sentinels), which the scan and the traceback read. After the encoding the letters
//! are read only by the batch's query chunks (`search_query_chunks`: their cut, DUST and encoding)
//! and by the reports that show a query's letters. With the switch, where the rule below holds,
//! `search_cli` moves the records into a `Store` and the batch loop takes each batch's records out
//! of it (`x_query_stream::Queries::Taken`), so the loop holds them alone: `search_query_batch`
//! drops their letters after its last reader (the encoding and the query chunks), before it builds
//! the lookup tables (`record.sequence = Vec::new()`, as LOSAT_X_BNSTORE does with the subjects'
//! letters). The loop keeps the records without letters for a writer after the batch, as the
//! streamed records are kept. A batch split into query chunks cuts them from its records' letters,
//! so it drops the letters after its chunks, before its traceback; the chunks' own copies
//! (`FastaRecord::cut`) keep theirs until the chunk ends, as before (dropping a copy of about 1 MiB
//! before its chunk's scan saved nothing measurable: blastn of the 5.5 Mb query below, whose
//! records LOSAT_X_QSTREAM reads batch by batch, read the same max RSS).
//!
//! The rule (`Rule::acts`; the reach of the switch is recorded where `search_cli` decides): every
//! output format reads no query letter after the search (`BlastnOutputFormat::reads_query_letters`:
//! the tabular formats 6 and 7, the Seq-aligns of 8, 9 and 12, and SAM without both `SQ` and `SR`;
//! not outfmt 0, 5, 11, 13 to 16, nor SAM's `SQ SR`) and no caller takes the hits' records (their
//! query rows), without `-query_loc` (its input records keep their letters), and for the records
//! read before the search only (not a streamed search, and not `run_local`, whose caller holds its
//! records).
//!
//! Why the reports are the same. The records keep their local IDs, titles and messages; every
//! length the search and the reports read after the encoding comes from the lengths taken before
//! it (`query_lengths` in the traceback, the loop's searched lengths for an invalid query's
//! message, `QueryNames::searched_lengths` in the Seq-aligns), and no format of the rule reads a
//! query letter (the HSP rows of outfmt 0 and the XML and JSON formats, the archive's Bioseqs and
//! its DUST, and SAM's SEQ of the query with `SQ SR` do). Dropping memory changes no value that the
//! search computes. No floating-point operation is involved.
//!
//! Memory: the peak of a megablast batch is the scan (with BNLUTDROP, graduated), where the letters
//! of a 5.5 Mb query (E. coli Sakai x MG1655) were 5.2 MiB beside the 10.5 MiB of the BLASTNA
//! buffer and the 106 MiB of the lookup table (plan v4 line MF, `REPORT.md` §2.3): its peak
//! `RssAnon` falls by 4.4-6.4 MiB at one and 16 threads, `-subject` and `-db`. Megablast of 31,168
//! cDNAs at one thread (37 batches, 89 MB of letters) drops each batch's letters after its
//! encoding; its peak, in the first batches, falls by 4.4 MiB.
//!
//! Timing: with LOSAT_TIMING, `[TIMING] x_bn_q_drop: batches= records= letters=` gives the number
//! of batches and records whose letters were dropped, and their letters.
//!
//! The letters packed. Dropping a batch's letters after its encoding alone frees them late: the
//! letters of every record were touched while the file was read, and the pages stay with the process
//! (the heap reuses them for later allocations, so the peak falls by little; plan v4 line A10 session
//! 19, the `RssAnon` trace of 31,168 cDNAs). So where the store acts (the rule above), the reader
//! packs each nucleotide record as it assembles it (`FastaInputSource::x_pack_letters`: 2-bit codes,
//! the runs of other bytes, the lowercase mask as intervals, `sequence::packed_letters`) in place of
//! the copy of its buffer, which the next record reuses; the store keeps the packed letters beside the records
//! (`Store::new`), and the loop makes a batch's letters again when it takes the batch (`Store::take`,
//! by the search pool's threads with two or more), before any reader of them (WindowMasker, DUST, the
//! encoding). The letters made again are the reader's, byte for byte (`U` read as `T`, the lowercase
//! mask in lower case); a record whose store would hold no fewer bytes than its letters keeps them.
//! No search value changes. With LOSAT_TIMING, `[TIMING] x_bn_q_pack: records= letters= bytes=
//! unpack_s=` gives the records packed, their letters, the bytes held for them and the time spent
//! making their letters again.
//!
//! `LOSAT_X_BNQSLOTS`, which opened the store to the query-batch slots of LOSAT_X_QBATCHPAR (the
//! store lent the loop's batch to a slot and copied the records of batches ahead), was removed with
//! the slots on 2026-10-28.
//!
//! The store in parts (EXPERIMENT LOSAT_X_BNQPARTS; in the strict set from the pull request that
//! added it, on unless it is `0`; it acts wherever the store does). The store above keeps each
//! record the reader made (its local ID, title and packed letters on the heap, beside the 104
//! bytes of the record), all read on the calling thread before the pool exists and all freed on one
//! thread at the end. With the switch the store holds parts (`Part`): one buffer per part with each
//! record's title, messages and letters (their store, `PackedLetters::append_to`, where the reader
//! packed them; else the letters) and a table of where each record ends. The reader's record is
//! dropped once it is copied, so the next record reuses its allocations. The parts are read on the
//! search's threads where the file can be cut at record starts as LOSAT_X_QSTREAM's first reading
//! is cut (`x_query_scan_par::store_parts`), else by the search's reader into parts of about 1 MiB
//! (`read_parts`). When the loop takes a batch, its records are made again from their parts
//! (`Part::record`, by the pool's threads with two or more): the title and messages copied, the
//! letters made again (`unpack_from`), the local ID `Query_<n>` from the record's position among all
//! records, which is the ID the reader's counter gave it (`Part::push` keeps any other ID, the ones
//! `-parse_deflines` read, beside the part with the record's source). A part is dropped once the
//! loop has taken all its records. No search value changes. With `LOSAT_X_STATS`,
//! `[X_STATS] BNQPARTS records= parts= in_parts=` (in_parts: read on the search's threads).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

use crate::algorithm::blastn::coordination::QueryLetters;
use crate::algorithm::blastn::hsp::BlastnOutputFormat;
use crate::blastinput::fasta_reader::{
    FastaInputSource, FastaRecord, ParseErrorCode, QueryEnd, ReadError, RecordSource,
};
use crate::blastinput::seq_range::QueryInput;
use crate::sequence::packed_letters::{unpack_from, PackedAt, PackedLetters};
use crate::utils::threading::SearchPool;

/// No NCBI counterpart: LOSAT_X_BNQPACK (read once).
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNQPACK"))
}

/// No NCBI counterpart: LOSAT_X_BNQPARTS (read once).
pub(super) fn parts_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNQPARTS"))
}

/// No NCBI counterpart: what `search_cli` knows of a search when it decides whether the switch
/// acts.
pub(super) struct Rule<'a> {
    pub formats: &'a [BlastnOutputFormat],
    /// `-query_loc` is given.
    pub query_range: bool,
    /// The caller takes the hits' records (`ReportOutputs::hits`, which hold the query rows).
    pub hit_records: bool,
}

impl Rule<'_> {
    /// Whether the records read before the search are moved into a `Store` (the switch on).
    pub(super) fn acts(&self) -> bool {
        !self.query_range
            && !self.hit_records
            && !self
                .formats
                .iter()
                .any(|format| format.reads_query_letters())
    }
}

/// The query records read before a search, which the batch loop takes out batch by batch, in
/// order (`x_query_stream::Queries::Taken`).
pub(super) struct Store {
    /// The records the loop has not taken, the packed letters of those packed, and how many
    /// records the loop has taken.
    left: Mutex<Left>,
    /// EXPERIMENT (LOSAT_X_BNQPARTS): the records held in parts (`Store::from_parts`; `left` is then
    /// empty).
    parts: Option<Mutex<PartsLeft>>,
    /// Every record searched whole (no `-query_loc`).
    query_input: QueryInput,
}

/// EXPERIMENT (LOSAT_X_BNQPARTS): the parts with records the loop has not taken, each with the
/// position of its first record among all records, and how many records the loop has taken.
struct PartsLeft {
    parts: std::collections::VecDeque<(usize, Part)>,
    taken: usize,
    /// The prefix of the local IDs the reader's counter gives (`ReaderConfig::id_prefix`).
    id_prefix: &'static str,
    /// The parts and whether they were read on the search's threads, for the `[X_STATS]` line.
    count: usize,
    in_parts: bool,
}

struct Left {
    records: std::vec::IntoIter<FastaRecord>,
    /// The packed letters of the records not taken, with their positions among all records.
    packed: std::collections::VecDeque<(usize, Arc<PackedLetters>)>,
    taken: usize,
}

/// The letters of a batch's records for the readers that set the batch up
/// (`coordination::QueryLetters`: WindowMasker, DUST, the lower-case masks, the encoding, the query
/// chunks' cut): each record's own letters.
#[derive(Clone, Copy)]
pub(super) struct Letters<'a> {
    records: &'a [FastaRecord],
}

impl Letters<'_> {
    /// Record `index` cut to `from..to`, with its letters (a query chunk, `FastaRecord::cut`).
    pub(super) fn cut(&self, index: usize, from: usize, to: usize) -> FastaRecord {
        self.records[index].cut(from, to)
    }
}

impl QueryLetters for Letters<'_> {
    fn count(&self) -> usize {
        self.records.len()
    }

    fn length(&self, index: usize) -> usize {
        self.records[index].sequence.len()
    }

    fn read<T>(&self, index: usize, read: impl FnOnce(&[u8]) -> T) -> T {
        read(&self.records[index].sequence)
    }
}

impl Store {
    /// The store of `records`, and of the letters of those packed as they were read (`packed`:
    /// their positions among the records, in order, and their letters; empty when none was).
    pub(super) fn new(records: Vec<FastaRecord>, packed: Vec<(usize, PackedLetters)>) -> Self {
        let mut query_input = QueryInput::whole(&records);
        for (index, store) in &packed {
            query_input.input_lengths[*index] = store.len();
        }
        if !packed.is_empty() {
            PACKED.fetch_add(packed.len() as u64, Ordering::Relaxed);
            PACKED_LETTERS.fetch_add(
                packed.iter().map(|(_, store)| store.len() as u64).sum(),
                Ordering::Relaxed,
            );
            PACKED_BYTES.fetch_add(
                packed
                    .iter()
                    .map(|(_, store)| store.heap_bytes() as u64)
                    .sum(),
                Ordering::Relaxed,
            );
        }
        Self {
            query_input,
            left: Mutex::new(Left {
                records: records.into_iter(),
                packed: packed
                    .into_iter()
                    .map(|(index, store)| (index, Arc::new(store)))
                    .collect(),
                taken: 0,
            }),
            parts: None,
        }
    }

    /// EXPERIMENT (LOSAT_X_BNQPARTS): the store of the records of `parts`, in order (`in_parts`:
    /// read on the search's threads, `x_query_scan_par::store_parts`, else by the search's reader);
    /// `id_prefix` is the reader's (`ReaderConfig::id_prefix`).
    pub(super) fn from_parts(parts: Vec<Part>, id_prefix: &'static str, in_parts: bool) -> Self {
        let count = parts.iter().map(Part::len).sum();
        let mut input_lengths = Vec::with_capacity(count);
        let (mut packed, mut packed_letters, mut bytes) = (0, 0, 0);
        for part in &parts {
            for record in &part.records {
                input_lengths.push(record.letters.len());
                if let PartLetters::Packed(at) = record.letters {
                    packed += 1;
                    packed_letters += at.len() as u64;
                    bytes += at.bytes() as u64;
                }
            }
        }
        if packed > 0 {
            PACKED.fetch_add(packed, Ordering::Relaxed);
            PACKED_LETTERS.fetch_add(packed_letters, Ordering::Relaxed);
            PACKED_BYTES.fetch_add(bytes, Ordering::Relaxed);
        }
        let mut first = 0;
        let parts_count = parts.len();
        let parts = parts
            .into_iter()
            .map(|part| {
                let at = first;
                first += part.len();
                (at, part)
            })
            .collect();
        Self {
            query_input: QueryInput {
                placements: Default::default(),
                ordinals: (0..count).collect(),
                input_lengths,
                skipped: vec![false; count],
                range_length: None,
            },
            left: Mutex::new(Left {
                records: Vec::new().into_iter(),
                packed: Default::default(),
                taken: 0,
            }),
            parts: Some(Mutex::new(PartsLeft {
                parts,
                taken: 0,
                id_prefix,
                count: parts_count,
                in_parts,
            })),
        }
    }

    /// The query input (every record searched whole).
    pub(super) fn query_input(&self) -> &QueryInput {
        &self.query_input
    }

    /// The length of every record.
    pub(super) fn lengths(&self) -> &[usize] {
        &self.query_input.input_lengths
    }

    /// The records `range`, taken out: the loop takes its batches in order, each starting where
    /// the one before ended. Packed letters are made again (by `pool`'s threads with two or more).
    pub(super) fn take(
        &self,
        range: std::ops::Range<usize>,
        pool: &SearchPool<'_>,
    ) -> anyhow::Result<Vec<FastaRecord>> {
        if let Some(parts) = &self.parts {
            let mut parts = parts.lock().unwrap_or_else(PoisonError::into_inner);
            return parts.take(range, self.lengths().len(), pool);
        }
        let mut left = self.left.lock().unwrap_or_else(PoisonError::into_inner);
        let Left {
            records,
            packed,
            taken,
        } = &mut *left;
        anyhow::ensure!(
            range.start == *taken && range.len() <= records.len(),
            "LOSAT_X_BNQPACK: the batch of the records {range:?} after those up to {taken}"
        );
        *taken = range.end;
        let mut batch: Vec<FastaRecord> = records.by_ref().take(range.len()).collect();
        let mut stores: Vec<Option<Arc<PackedLetters>>> = Vec::new();
        while packed.front().is_some_and(|(index, _)| *index < range.end) {
            let (index, store) = packed.pop_front().expect("a front");
            stores.resize_with(range.len(), || None);
            stores[index - range.start] = Some(store);
        }
        drop(left);
        unpack(&mut batch, stores, pool);
        Ok(batch)
    }

    /// Writes the reader's messages of every record (none taken yet).
    pub(super) fn write_warnings(
        &self,
        out: &mut (dyn std::io::Write + Send),
    ) -> anyhow::Result<()> {
        if let Some(parts) = &self.parts {
            let parts = parts.lock().unwrap_or_else(PoisonError::into_inner);
            for (first, part) in &parts.parts {
                for index in parts.taken.saturating_sub(*first)..part.len() {
                    out.write_all(part.warnings(index))?;
                }
            }
            return Ok(());
        }
        let left = self.left.lock().unwrap_or_else(PoisonError::into_inner);
        for record in left.records.as_slice() {
            out.write_all(&record.warnings)?;
        }
        Ok(())
    }
}

// No NCBI counterpart: the summary line of LOSAT_X_BNQPARTS (with LOSAT_X_STATS); it prints
// counters and changes no value NCBI computes.
impl Drop for Store {
    fn drop(&mut self) {
        let Some(parts) = &self.parts else {
            return;
        };
        if std::thread::panicking() || std::env::var_os("LOSAT_X_STATS").is_none() {
            return;
        }
        let parts = parts.lock().unwrap_or_else(PoisonError::into_inner);
        eprintln!(
            "[X_STATS] BNQPARTS records={} parts={} in_parts={}",
            self.query_input.input_lengths.len(),
            parts.count,
            u8::from(parts.in_parts)
        );
    }
}

impl PartsLeft {
    /// `Store::take` of the parts: the records `range` made again from their parts (by `pool`'s
    /// threads with two or more), and the parts all of whose records are taken dropped.
    fn take(
        &mut self,
        range: std::ops::Range<usize>,
        count: usize,
        pool: &SearchPool<'_>,
    ) -> anyhow::Result<Vec<FastaRecord>> {
        anyhow::ensure!(
            range.start == self.taken && range.end <= count,
            "LOSAT_X_BNQPARTS: the batch of the records {range:?} after those up to {}",
            self.taken
        );
        let t0 = crate::utils::xtiming::clock();
        // The runs of the batch's records in each part: (part, first record in the part, first
        // record in the batch).
        let mut runs = Vec::new();
        let mut at = range.start;
        for (index, (first, part)) in self.parts.iter().enumerate() {
            if at >= range.end {
                break;
            }
            if at < first + part.len() {
                runs.push((index, at - first, at - range.start));
                at = (first + part.len()).min(range.end);
            }
        }
        // The letters of the batch's records that have their own vector, moved out of their parts.
        let mut owned: Vec<Option<Vec<u8>>> = Vec::new();
        for (run, &(index, from, offset)) in runs.iter().enumerate() {
            let end = runs.get(run + 1).map_or(range.len(), |&(_, _, next)| next) - offset + from;
            for (at, letters) in &mut self.parts[index].1.own {
                if (from..end).contains(at) {
                    owned.resize_with(range.len(), || None);
                    owned[offset + *at - from] = Some(std::mem::take(letters));
                }
            }
        }
        let parts = &self.parts;
        let id_prefix = self.id_prefix;
        let make = |offset: usize, own: Option<Vec<u8>>| {
            let run = runs[runs.partition_point(|&(_, _, from)| from <= offset) - 1];
            let (first, part) = &parts[run.0];
            let index = run.1 + offset - run.2;
            part.record(index, first + index, id_prefix, own)
        };
        let batch: Vec<FastaRecord>;
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        {
            use rayon::prelude::*;
            batch = match (pool.enabled() && pool.threads() > 1, owned.is_empty()) {
                (true, true) => pool.install(|| {
                    (0..range.len())
                        .into_par_iter()
                        .map(|offset| make(offset, None))
                        .collect()
                }),
                (true, false) => pool.install(|| {
                    owned
                        .into_par_iter()
                        .enumerate()
                        .map(|(offset, own)| make(offset, own))
                        .collect()
                }),
                (false, true) => (0..range.len()).map(|offset| make(offset, None)).collect(),
                (false, false) => owned
                    .into_iter()
                    .enumerate()
                    .map(|(offset, own)| make(offset, own))
                    .collect(),
            };
        }
        #[cfg(not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )))]
        {
            let _ = pool;
            owned.resize_with(range.len(), || None);
            batch = owned
                .into_iter()
                .enumerate()
                .map(|(offset, own)| make(offset, own))
                .collect();
        }
        self.taken = range.end;
        while self
            .parts
            .front()
            .is_some_and(|(first, part)| first + part.len() <= self.taken)
        {
            self.parts.pop_front();
        }
        add_unpack_time(t0);
        Ok(batch)
    }
}

/// EXPERIMENT (LOSAT_X_BNQPARTS): the records of a run of consecutive records of the query file (a
/// part of it read on a thread of its own, `x_query_scan_par::store_parts`, or the whole file read
/// by the search's reader, `read_part`) held in one buffer: each record's title, reader messages
/// and letters (their store where the reader packed them, `PackedLetters::append_to`, else the
/// letters), made again when the loop takes their batch (`Part::record`).
#[derive(Default)]
pub(super) struct Part {
    bytes: Vec<u8>,
    records: Vec<PartRecord>,
    /// The records whose local ID is not the one their reader's counter gave (`-parse_deflines`)
    /// or that carry a source (`FastaRecord::source`): their position in the part, ID and source.
    special: Vec<(usize, String, Option<Box<RecordSource>>)>,
    /// The letters of the records that the reader did not pack and that have `OWN_LETTERS` or
    /// more, in their own vector (moved in and out, not copied), by their position in the part.
    own: Vec<(usize, Vec<u8>)>,
}

/// The letters of a record the reader did not pack from which the record keeps its own vector.
const OWN_LETTERS: usize = 1 << 16;

#[derive(Clone, Copy)]
struct PartRecord {
    /// Where the record ends in `bytes` (it starts where the one before ends).
    end: usize,
    title: usize,
    warnings: usize,
    letters: PartLetters,
}

#[derive(Clone, Copy)]
enum PartLetters {
    /// The letters themselves (the reader did not pack them).
    Letters(usize),
    /// The letters in their own vector (`Part::own`).
    Own(usize),
    Packed(PackedAt),
}

impl PartLetters {
    fn len(&self) -> usize {
        match self {
            PartLetters::Letters(length) | PartLetters::Own(length) => *length,
            PartLetters::Packed(at) => at.len(),
        }
    }
}

impl Part {
    /// An empty part with room for `bytes` bytes and `records` records.
    pub(super) fn with_capacity(bytes: usize, records: usize) -> Part {
        Part {
            bytes: Vec::with_capacity(bytes),
            records: Vec::with_capacity(records),
            special: Vec::new(),
            own: Vec::new(),
        }
    }

    /// The bytes of the buffer that `push` would take for `record` (with the reader's `packed`).
    pub(super) fn bytes_for(record: &FastaRecord, packed: Option<&PackedLetters>) -> usize {
        record.title.len()
            + record.warnings.len()
            + match packed {
                Some(store) => store.stored_bytes(),
                None if record.sequence.len() >= OWN_LETTERS => 0,
                None => record.sequence.len(),
            }
    }

    /// The bytes held in the buffer.
    pub(super) fn bytes(&self) -> usize {
        self.bytes.len()
    }

    /// Gives back the room no record took (a part is complete).
    pub(super) fn shrink(&mut self) {
        self.bytes.shrink_to_fit();
        self.records.shrink_to_fit();
    }

    /// The number of records.
    pub(super) fn len(&self) -> usize {
        self.records.len()
    }

    /// Adds `record`, the `n`th record its reader read (from 1), with the store of its letters
    /// where the reader packed them (`packed`; its `sequence` is then empty). Returns whether its
    /// local ID is the one the reader's counter gave (`id_prefix` and `n`) and it has no source.
    pub(super) fn push(
        &mut self,
        record: FastaRecord,
        packed: Option<PackedLetters>,
        id_prefix: &str,
        n: usize,
    ) -> bool {
        let FastaRecord {
            local_id,
            title,
            sequence,
            warnings,
            source,
        } = record;
        self.bytes.extend_from_slice(&title);
        self.bytes.extend_from_slice(&warnings);
        let letters = match packed {
            Some(store) => PartLetters::Packed(store.append_to(&mut self.bytes)),
            None if sequence.len() >= OWN_LETTERS => {
                let length = sequence.len();
                self.own.push((self.records.len(), sequence));
                PartLetters::Own(length)
            }
            None => {
                self.bytes.extend_from_slice(&sequence);
                PartLetters::Letters(sequence.len())
            }
        };
        let counted = source.is_none() && counted_id(&local_id, id_prefix, n);
        if !counted {
            self.special.push((self.records.len(), local_id, source));
        }
        self.records.push(PartRecord {
            end: self.bytes.len(),
            title: title.len(),
            warnings: warnings.len(),
            letters,
        });
        counted
    }

    /// The bytes of record `index`: its title, messages and letters (or their store).
    fn parts_of(&self, index: usize) -> (&[u8], &[u8], &[u8]) {
        let start = index
            .checked_sub(1)
            .map_or(0, |before| self.records[before].end);
        let record = &self.records[index];
        let (title, rest) = self.bytes[start..record.end].split_at(record.title);
        let (warnings, letters) = rest.split_at(record.warnings);
        (title, warnings, letters)
    }

    /// The reader's messages of record `index`.
    fn warnings(&self, index: usize) -> &[u8] {
        self.parts_of(index).1
    }

    /// Record `index` made again, the input record `ordinal` (from 0) of a reader whose local IDs
    /// start with `id_prefix`; `own` is its letters where they have their own vector (moved out of
    /// the part by the caller).
    fn record(
        &self,
        index: usize,
        ordinal: usize,
        id_prefix: &str,
        own: Option<Vec<u8>>,
    ) -> FastaRecord {
        let (title, warnings, letters) = self.parts_of(index);
        let (local_id, source) = match self.special.binary_search_by_key(&index, |(at, _, _)| *at) {
            Ok(found) => {
                let (_, local_id, source) = &self.special[found];
                (local_id.clone(), source.clone())
            }
            Err(_) => (format!("{id_prefix}{}", ordinal + 1), None),
        };
        let mut sequence = Vec::new();
        match self.records[index].letters {
            PartLetters::Letters(_) => sequence.extend_from_slice(letters),
            PartLetters::Own(_) => {
                sequence = own.expect("the loop takes a record's own letters once");
            }
            PartLetters::Packed(at) => unpack_from(letters, at, &mut sequence),
        }
        FastaRecord {
            local_id,
            title: title.to_vec(),
            sequence,
            warnings: warnings.to_vec(),
            source,
        }
    }
}

/// Whether `id` is the local ID that a reader's counter gives its `n`th record (`generate_id`:
/// the prefix and `n` in decimal).
fn counted_id(id: &str, id_prefix: &str, n: usize) -> bool {
    id.strip_prefix(id_prefix).is_some_and(|digits| {
        !digits.starts_with('0')
            && digits.bytes().all(|byte| byte.is_ascii_digit())
            && digits.parse::<usize>() == Ok(n)
    })
}

/// EXPERIMENT (LOSAT_X_BNQPARTS): `x_read_queries_packed` into parts of at most `PART_BYTES` bytes
/// (a record that takes more has a part of its own): the records of `source` (their letters packed
/// as they are read) and how the reading ended, as `read_queries` reads them. A part is closed before
/// a record that would pass the mark, so its buffer is not copied as it grows, and then gives back
/// the room it did not take.
pub(super) fn read_parts<R: std::io::Read>(
    source: &mut FastaInputSource<R>,
    id_prefix: &str,
) -> (Vec<Part>, QueryEnd) {
    read_parts_of(source, id_prefix, PART_BYTES)
}

/// `read_parts` with parts closed at `part_bytes` (the tests make them small).
fn read_parts_of<R: std::io::Read>(
    source: &mut FastaInputSource<R>,
    id_prefix: &str,
    part_bytes: usize,
) -> (Vec<Part>, QueryEnd) {
    source.x_pack_letters();
    let mut parts = Vec::new();
    let mut part = Part::with_capacity(part_bytes, 0);
    let mut read = 0;
    let end = loop {
        if source.end() {
            break QueryEnd::Input;
        }
        let mut warnings = Vec::new();
        match source.next_sequence(&mut |message: &[u8]| {
            warnings.extend_from_slice(message);
            Ok(())
        }) {
            Ok(record) => {
                let packed = source.x_take_packed();
                read += 1;
                if part.len() > 0
                    && part.bytes() + Part::bytes_for(&record, packed.as_ref()) > part_bytes
                {
                    let mut done = std::mem::replace(&mut part, Part::with_capacity(part_bytes, 0));
                    done.shrink();
                    parts.push(done);
                }
                part.push(record, packed, id_prefix, read);
            }
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => break QueryEnd::EofError,
            Err(error) => break QueryEnd::Error { error, warnings },
        }
    };
    if part.len() > 0 || parts.is_empty() {
        part.shrink();
        parts.push(part);
    }
    (parts, end)
}

/// The most bytes of a part of `read_parts` (but for a record that takes more).
const PART_BYTES: usize = 1 << 20;

/// No NCBI counterpart: the letters of a batch's packed records made again, by the search pool
/// where it has two threads or more (each record by one thread).
fn unpack(
    batch: &mut [FastaRecord],
    stores: Vec<Option<Arc<PackedLetters>>>,
    pool: &SearchPool<'_>,
) {
    if stores.is_empty() {
        return;
    }
    let t0 = crate::utils::xtiming::clock();
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() && pool.threads() > 1 {
        use rayon::prelude::*;
        pool.install(|| {
            batch
                .par_iter_mut()
                .zip(stores.into_par_iter())
                .for_each(|(record, store)| {
                    if let Some(store) = store {
                        store.unpack_into(&mut record.sequence);
                    }
                })
        });
        add_unpack_time(t0);
        return;
    }
    let _ = pool;
    for (record, store) in batch.iter_mut().zip(stores) {
        if let Some(store) = store {
            store.unpack_into(&mut record.sequence);
        }
    }
    add_unpack_time(t0);
}

fn add_unpack_time(t0: Option<std::time::Instant>) {
    if let Some(t0) = t0 {
        UNPACK_NANOS.fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }
}

/// The query records of a batch as `search_query_batch` reads them: shared with other readers, or
/// held by the caller alone (the loop's batch taken from a `Store`, whose letters were made when it
/// was taken), whose letters the batch drops after its last reader of them.
pub(super) enum BatchQueries<'r> {
    Shared(&'r [FastaRecord]),
    Held(&'r mut [FastaRecord]),
}

impl BatchQueries<'_> {
    pub(super) fn records(&self) -> &[FastaRecord] {
        match self {
            BatchQueries::Shared(records) => records,
            BatchQueries::Held(records) => records,
        }
    }

    /// The records' letters for the readers that set the batch up (`Letters`).
    pub(super) fn letters(&self) -> Letters<'_> {
        match self {
            BatchQueries::Shared(records) => Letters { records },
            BatchQueries::Held(records) => Letters { records },
        }
    }

    /// The batch has encoded its letters and searched its query chunks: the letters of records the
    /// caller holds alone are dropped (with the switch).
    pub(super) fn after_encoding(&mut self) {
        let BatchQueries::Held(records) = self else {
            return;
        };
        if !on() {
            return;
        }
        let mut letters = 0;
        for record in records.iter_mut() {
            letters += record.sequence.len() as u64;
            record.sequence = Vec::new();
        }
        BATCHES.fetch_add(1, Ordering::Relaxed);
        RECORDS.fetch_add(records.len() as u64, Ordering::Relaxed);
        LETTERS.fetch_add(letters, Ordering::Relaxed);
    }
}

// No NCBI counterpart: the batches, records and letters dropped after their encoding.
static BATCHES: AtomicU64 = AtomicU64::new(0);
static RECORDS: AtomicU64 = AtomicU64::new(0);
static LETTERS: AtomicU64 = AtomicU64::new(0);
// No NCBI counterpart: the records packed as they were read, their letters and
// the bytes held for them, and the time spent making their letters again (with LOSAT_TIMING).
static PACKED: AtomicU64 = AtomicU64::new(0);
static PACKED_LETTERS: AtomicU64 = AtomicU64::new(0);
static PACKED_BYTES: AtomicU64 = AtomicU64::new(0);
static UNPACK_NANOS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: with LOSAT_TIMING, the `[TIMING] x_bn_q_drop` line `main` prints at exit
/// (when a batch dropped its letters).
pub(super) fn print_summary() {
    let batches = BATCHES.load(Ordering::Relaxed);
    if crate::utils::xtiming::on() && batches > 0 {
        eprintln!(
            "[TIMING] x_bn_q_drop: batches={batches} records={} letters={}",
            RECORDS.load(Ordering::Relaxed),
            LETTERS.load(Ordering::Relaxed)
        );
    }
    let packed = PACKED.load(Ordering::Relaxed);
    if crate::utils::xtiming::on() && packed > 0 {
        eprintln!(
            "[TIMING] x_bn_q_pack: records={packed} letters={} bytes={} unpack_s={:.3}",
            PACKED_LETTERS.load(Ordering::Relaxed),
            PACKED_BYTES.load(Ordering::Relaxed),
            UNPACK_NANOS.load(Ordering::Relaxed) as f64 / 1e9
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::fasta_reader::{read_queries, ReaderConfig};

    /// Nine records (three batches of three), every second one packed, as `search_cli` stores them.
    fn store() -> (Vec<FastaRecord>, Store) {
        let records: Vec<FastaRecord> = (0..9)
            .map(|index| {
                // Mostly ACGT (packed), a run of N and a lower-case stretch (kept beside the
                // codes).
                let mut letters: Vec<u8> = b"ACGTTGCAAGCTTCGAGATCCTAGGCATTACG"
                    .iter()
                    .cycle()
                    .skip(index)
                    .take(400 + 37 * index)
                    .copied()
                    .collect();
                letters[100..104].copy_from_slice(b"NNNN");
                letters[200..240].make_ascii_lowercase();
                FastaRecord::new(format!("Query_{index}"), b"title", &letters)
            })
            .collect();
        let mut held = records.clone();
        let mut packed = Vec::new();
        for (index, record) in held
            .iter_mut()
            .enumerate()
            .filter(|(index, _)| index % 2 == 1)
        {
            if let Some(store) = PackedLetters::pack(&record.sequence) {
                record.sequence = Vec::new();
                packed.push((index, store));
            }
        }
        assert!(!packed.is_empty(), "no record packed");
        (records, Store::new(held, packed))
    }

    /// The records of `text` as `read_queries` reads them, and how the reading ended.
    fn read_all(text: &[u8], config: ReaderConfig) -> (Vec<FastaRecord>, String) {
        let records = read_queries(&mut FastaInputSource::from_bytes(text, config));
        (records.records, format!("{:?}", records.end))
    }

    /// The records of `text` read into parts of `part_bytes` and taken back out of their store in
    /// batches of `batch` records (by a pool of `threads` threads), how the reading ended and the
    /// number of parts.
    fn through_parts(
        text: &[u8],
        config: ReaderConfig,
        part_bytes: usize,
        batch: usize,
        threads: usize,
    ) -> (Vec<FastaRecord>, String, usize) {
        let mut source = FastaInputSource::from_bytes(text, config);
        let (parts, end) = read_parts_of(&mut source, config.id_prefix, part_bytes);
        let count = parts.len();
        let store = Store::from_parts(parts, config.id_prefix, false);
        let total = store.lengths().len();
        let take = |pool: &SearchPool<'_>| -> anyhow::Result<Vec<FastaRecord>> {
            let mut records = Vec::new();
            for from in (0..total).step_by(batch.max(1)) {
                records.extend(store.take(from..(from + batch).min(total), pool)?);
            }
            Ok(records)
        };
        let records = if threads > 1 {
            crate::utils::threading::with_search_pool(threads, "blastn", take).expect("take")
        } else {
            take(&SearchPool::x_serial()).expect("take")
        };
        (records, format!("{end:?}"), count)
    }

    /// FASTA of `records` records: titles, `N` runs, lower case, a record without letters, one
    /// without an id, ids of several forms (for `-parse_deflines`), then `tail`.
    fn text(records: usize, tail: &[u8]) -> Vec<u8> {
        let mut text = Vec::new();
        for index in 0..records {
            match index % 7 {
                0 => text.extend_from_slice(b">\n"),
                1 => {
                    text.extend_from_slice(format!(">gi|{index}|ref|NC_{index}.1| t\n").as_bytes())
                }
                2 => text.extend_from_slice(format!(">lcl|q{index}\n").as_bytes()),
                _ => text.extend_from_slice(format!(">read{index} a title {index}\n").as_bytes()),
            }
            if index % 11 == 5 {
                continue;
            }
            let letters: Vec<u8> = b"ACGTTGCAAGCTTCGAGATCCTAGGCATTACG"
                .iter()
                .cycle()
                .skip(index)
                .take(20 + (index * 37) % 300)
                .copied()
                .collect();
            for line in letters.chunks(60) {
                text.extend_from_slice(line);
                text.push(b'\n');
            }
            if index % 9 == 4 {
                text.extend_from_slice(b"NNNNNNNNNNacgtacgtacgtnnnnACGT\n");
            }
        }
        text.extend_from_slice(tail);
        text
    }

    // The records read into parts and taken back out of their store are the reader's, with how the
    // reading ended (the end, a blank-line `eEOF`, a reader error) and the reader's messages, also
    // with the ids `-parse_deflines` reads, across parts and batches of any size, by one thread and
    // the pool.
    #[test]
    fn records_read_into_parts_are_the_readers_records() {
        let plain = ReaderConfig::query("BLASTN", false, false);
        let parsed = plain.with_parse_deflines(true);
        // A record that the reader does not pack (a run on almost every letter) and that is long
        // enough to keep its own vector.
        let mut dense = b">dense runs\n".to_vec();
        for line in b"ARCYGKTMNSAW".repeat(6_000).chunks(70) {
            dense.extend_from_slice(line);
            dense.push(b'\n');
        }
        dense.extend_from_slice(b">after\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        for (tail, ends) in [
            (dense.as_slice(), "Input"),
            (b"".as_slice(), "Input"),
            (b"\n; a comment\n\n".as_slice(), ""),
            (b">hyphen\nACGT-ACGT\n".as_slice(), "Input"),
            (b">bad\n%%%%%%%%%%%%%%%%\n".as_slice(), "Error"),
        ] {
            for config in [plain, parsed] {
                let text = text(400, tail);
                let (expected, end) = read_all(&text, config);
                assert!(end.starts_with(ends), "{end}");
                for part_bytes in [1, 700, 5_000, 1 << 20] {
                    for (batch, threads) in [(1, 1), (37, 1), (64, 4), (500, 3)] {
                        let (records, got_end, parts) =
                            through_parts(&text, config, part_bytes, batch, threads);
                        assert_eq!(got_end, end);
                        assert_eq!(records, expected, "part_bytes {part_bytes} batch {batch}");
                        if tail == dense.as_slice() {
                            let long = &records[records.len() - 2];
                            assert!(long.sequence.len() >= OWN_LETTERS);
                            assert!(PackedLetters::pack(&long.sequence).is_none());
                        }
                        assert!(part_bytes >= 5_000 || parts > 5, "{parts} parts");
                    }
                }
            }
        }
    }

    // The loop's batches taken out of the store are the records read, the packed letters made
    // again (one thread and the pool).
    #[test]
    fn the_loop_takes_its_batches_with_their_letters() {
        let (records, held) = store();
        let serial = SearchPool::x_serial();
        for range in [0..3, 3..6, 6..9] {
            let batch = held.take(range.clone(), &serial).expect("take");
            assert_eq!(&batch[..], &records[range]);
        }
        assert!(held.take(9..10, &serial).is_err());
    }
}
