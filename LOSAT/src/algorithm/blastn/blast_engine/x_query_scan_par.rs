//! EXPERIMENT (LOSAT_X_QSCANPAR, in the strict set: on unless it is `0`): the first reading of
//! LOSAT_X_QSTREAM (`x_query_stream::scan`, which keeps each query record's length and how the
//! reading ends) cut into parts at record starts and read on `-num_threads` threads (the calling
//! thread and scoped helpers, before the search pool exists), each part by a reader of its own;
//! and, where LOSAT_X_QSTREAM does not act, the reading of every record before the search
//! (`read_queries`, `x_read_queries_packed`) cut and read the same way (`records`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//! 	    BLAST_PROF_START( APP.LOOP.PRE );
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//! ```
//! NCBI reads its query batches one at a time, inside the loop. LOSAT needs every record's length
//! before its first batch (the batch planner, LOSAT_X_QBATCHPAR's plan, LOSAT_X_RPTSTREAM's checks),
//! so where LOSAT_X_QSTREAM acts it reads the whole query file once for the lengths alone, on one
//! thread, while every other thread of the search waits: for dc-megablast of 31,168 D. melanogaster
//! cDNAs (96 MB) against 3 Mb at 16 threads on the workstation, 0.12 s of a 1.55 s search. Where it
//! does not act (`-task megablast`, and the threaded WebAssembly command, whose WASI has no second
//! handle on a file), every record is read there, on one thread: 0.28-0.30 s of the same search's
//! 2.6 s in the threaded command (Node, workstation), before the first thread starts. No NCBI
//! counterpart: it is speed only.
//!
//! What it does. The bytes from the reader's start to the end of the file are cut at `parts`
//! targets (at most `MAX_SPAN` bytes apart, at least four per thread). Each part runs from the
//! first record start at or after its target (a `>` at the start of a line, not followed by `?`)
//! to the next part's. The threads take the parts in turn; each part is read by a new reader over
//! its bytes (`FastaInputSource::x_from_reader`, the same configuration) with the loop of `scan`:
//! `end()`, then `next_sequence`, the length of each record kept and the record dropped. The
//! lengths of the parts, in order, are the result. `records` keeps each record (and, with
//! LOSAT_X_BNQPACK, its packed letters, as `x_read_queries_packed` takes them) instead.
//!
//! Why the lengths are the first reading's. `scan` keeps only the records' lengths and how the
//! reading ends; the warnings of the records it reads are dropped there (the second reading posts
//! them per batch). A record's letters come from its lines alone: `read_fasta_seq` resets every
//! field of the record at its start, and what the reader carries from record to record is the
//! local-ID counter (`GenerateID`: ids only), the line number (messages and title lines only), the
//! end-of-line style and the pushback buffers of the line reader. A line that starts with `>` and
//! not with `>?` always starts a record: at the start of a record it is the defline, after the
//! defline it ends the record (`read_fasta_seq`'s second defline, ungot and read as the next
//! record's), and `>?_` is read as `>` only when what follows is not `?`. So at a part's first byte
//! the first reading is between two records, at the start of a line, with its next line that
//! defline, and a new reader started there reads the same lines into the same records, provided
//! the lines are cut where the first reading cuts them. They are where the file has no `\r`: then
//! every end-of-line style splits at `\n` alone (`Unknown` and `CrLf` look for `\r` too and never
//! find it), and the line reader pushes nothing back (its only pushback is a `\r` found inside a
//! line, `advance_simple`), so neither the buffered window of the input nor the bytes before the
//! part change a line. A part read to its end (its stream ends right after the line before the next
//! part's defline, so its last record ends there, as at that defline) gives the first reading's
//! records of those bytes. Only part 0 starts where the first reading starts, so only there can the
//! first line of the file be a Seq-id line (`read_one_seq` with the data loaders) or the reading
//! meet an `eEOF` before a record; a later part's first line is its defline.
//!
//! Why the records are those of one reading. Of what the reader carries, the records keep two:
//! the local ID (`Query_N`, the counter's N) and the messages, whose line numbers count from the
//! reader's first line. A part's reader starts its counter at 1, so `records` takes a part only
//! where each of its records has the ID of its place in the part (one count per record, as in one
//! reading of those lines) and gives it the ID of its place among all the records; and it takes a
//! part after the first only where none of its records has a message (the first part's messages
//! count lines as one reading does: it starts where the reading starts). Elsewhere every record is
//! read as before.
//!
//! Where the parts read. On unix each thread reads the search's handle at its offsets (`read_at`,
//! which leaves the file position alone). WASI preview 1 has neither such a read in Rust's stable
//! library nor a second handle (`try_clone`), so there each thread opens the path again and seeks
//! it, and takes no part unless the file it opened is a regular file of the handle's length and
//! modification time.
//!
//! When it does not act (the reading runs as before, on the search's reader, which no part has
//! moved): the switch is off; one thread; a target with neither unix file reads at an offset nor
//! WASI, or without the `parallel` feature (the serial WebAssembly command, Windows); standard
//! input or a file that is not a regular one; the BLAST archive (its title check,
//! `check_titles_for`, comes from the first reading); a file of fewer than `MIN_BYTES` bytes after
//! the reader's start (`LOSAT_X_QSCANPARMIN`, tuning, sets it); and, found while the parts are
//! read, a `\r` anywhere in the file, a `>?` anywhere in it (a gap line, whose length becomes that
//! many `N`s, or a `>?_` defline), a target with no record start within `MAX_REACH` bytes after it
//! (a record that long would take that much memory in every thread that reads one), a part whose
//! bytes cannot all be read (an I/O error or a file that ends early: a reader takes either for the
//! end of its input), a part whose reading ends with anything but the end of its bytes (an `eEOF`,
//! a reader error or LOSAT's rejection of a Seq-id line: the first reading then finds it, with its
//! messages, as before), and, for `records`, a part whose IDs or messages are not those of one
//! reading (above). What is read twice then is the parts' bytes up to the failure.
//!
//! Memory: for the lengths, a part holds at most `MAX_SPAN + MAX_REACH` bytes of records, and
//! without gap lines a record has no more letters than bytes, so a thread holds one record of at
//! most that many letters and the reader's buffers at a time. `records` holds what the reading
//! before the search holds (every record), and the parts' records until they are joined. No
//! floating-point operation is involved.
//!
//! The decision is marked with `x_reach!("LOSAT_X_QSCANPAR", acts)` where LOSAT_X_QSTREAM reads
//! the lengths, or the records are read before the search, and the switch is on, without the BLAST
//! archive (`acts`: the parts are taken).

#![cfg_attr(
    not(all(feature = "parallel", any(unix, target_os = "wasi"))),
    allow(dead_code)
)]

use std::fs::File;
use std::path::Path;
use std::sync::OnceLock;

use crate::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
use crate::sequence::packed_letters::PackedLetters;

/// The bytes after the reader's start below which the first reading stays on one thread: about
/// 4 ms of it on the workstation.
const MIN_BYTES: u64 = 4 << 20;
/// The most bytes between two targets.
const MAX_SPAN: u64 = 256 << 10;
/// The most bytes after a target searched for a record start.
const MAX_REACH: u64 = 256 << 10;
/// The parts per thread, at least.
const PARTS_PER_THREAD: u64 = 4;

/// Whether LOSAT_X_QSCANPAR is on (a strict-set switch: on unless it is `0`).
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_QSCANPAR"))
}

/// `LOSAT_X_QSCANPARMIN` (tuning): the bytes below which the first reading stays on one thread.
fn min_bytes() -> u64 {
    static MIN: OnceLock<u64> = OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_QSCANPARMIN")
            .ok()
            .and_then(|value| value.trim().parse::<u64>().ok())
            .unwrap_or(MIN_BYTES)
    })
}

/// The query file the parts read: `file`, a handle on it (the search's, or a second one), opened
/// at `path` (where the parts read, the module's header).
#[derive(Clone, Copy)]
pub(super) struct Input<'a> {
    pub(super) file: &'a File,
    pub(super) path: &'a Path,
}

/// No NCBI counterpart (NCBI reads the query once): a handle for `records` on the query file
/// `file` that the search opened at `path`, and the position its reader starts from, where `file`
/// is a regular file: a clone of `file` (`try_clone`, the same open file), or on WASI, which has
/// no descriptor duplication, the path opened again where it is a regular file of `file`'s length
/// and modification time.
pub(super) fn second_handle(path: &Path, file: &File) -> Option<(File, u64)> {
    use std::io::Seek;

    let metadata = file.metadata().ok()?;
    if !metadata.is_file() {
        return None;
    }
    let mut handle = file;
    let start = handle.stream_position().ok()?;
    let again = match file.try_clone() {
        Ok(again) => again,
        Err(_) if cfg!(target_os = "wasi") => {
            let again = File::open(path).ok()?;
            let reopened = again.metadata().ok()?;
            let same = reopened.is_file()
                && reopened.len() == metadata.len()
                && reopened.modified().ok()? == metadata.modified().ok()?;
            if !same {
                return None;
            }
            again
        }
        Err(_) => return None,
    };
    Some((again, start))
}

/// The lengths of the query records from byte `start` of `input` to its end, read in parts by
/// `threads` threads, or `None` where the first reading runs as before (the module's header).
/// `config` is the first reading's configuration, without the BLAST archive.
pub(super) fn lengths(
    input: Input<'_>,
    start: u64,
    config: ReaderConfig,
    threads: usize,
) -> Option<Vec<usize>> {
    debug_assert!(!config.keep_defline);
    let lengths = lengths_with(input, start, config, threads, min_bytes(), Limits::DEFAULT);
    crate::utils::xreach::x_reach!("LOSAT_X_QSCANPAR", lengths.is_some());
    lengths
}

/// The records of the query file from byte `start` of `input` to its end, as the reading before
/// the search reads them (`read_queries`, or with `pack` `x_read_queries_packed`: the records, a
/// packed one without letters, and the packed letters with their records' places), read in parts
/// by `threads` threads; `None` where that reading runs as before (the module's header), which
/// then also reads how the input ends. `config` is that reading's configuration, without the
/// BLAST archive. A reading of the parts ends at the end of the input (`QueryEnd::Input`).
pub(super) fn records(
    input: Input<'_>,
    start: u64,
    config: ReaderConfig,
    threads: usize,
    pack: bool,
) -> Option<PartRecords> {
    debug_assert!(!config.keep_defline);
    let records = records_with(
        input,
        start,
        config,
        threads,
        pack,
        min_bytes(),
        Limits::DEFAULT,
    );
    crate::utils::xreach::x_reach!("LOSAT_X_QSCANPAR", records.is_some());
    records
}

/// The records of `records`: every record in order, and with `pack` the packed letters of the
/// packed ones with their places among the records.
#[derive(Debug, Default)]
pub(super) struct PartRecords {
    pub(super) records: Vec<FastaRecord>,
    pub(super) packed: Vec<(usize, PackedLetters)>,
}

/// The sizes of the parts (tests make them small).
#[derive(Clone, Copy, Debug)]
struct Limits {
    max_span: u64,
    max_reach: u64,
    parts_per_thread: u64,
}

impl Limits {
    const DEFAULT: Limits = Limits {
        max_span: MAX_SPAN,
        max_reach: MAX_REACH,
        parts_per_thread: PARTS_PER_THREAD,
    };
}

fn lengths_with(
    input: Input<'_>,
    start: u64,
    config: ReaderConfig,
    threads: usize,
    min_bytes: u64,
    limits: Limits,
) -> Option<Vec<usize>> {
    let parts = in_parts(input, start, threads, min_bytes, limits, |at, part| {
        read_part(at, part, config, &mut |record| record.sequence.len())
    })?;
    Some(parts.into_iter().flatten().collect())
}

fn records_with(
    input: Input<'_>,
    start: u64,
    config: ReaderConfig,
    threads: usize,
    pack: bool,
    min_bytes: u64,
    limits: Limits,
) -> Option<PartRecords> {
    let parts = in_parts(input, start, threads, min_bytes, limits, |at, part| {
        read_part_records(at, part, config, pack)
    })?;
    let mut joined = PartRecords {
        records: Vec::with_capacity(parts.iter().map(|part| part.records.len()).sum()),
        packed: Vec::with_capacity(parts.iter().map(|part| part.packed.len()).sum()),
    };
    for (index, part) in parts.into_iter().enumerate() {
        let offset = joined.records.len();
        for mut record in part.records {
            if index > 0 {
                record.local_id = local_id(config, joined.records.len());
            }
            joined.records.push(record);
        }
        joined.packed.extend(
            part.packed
                .into_iter()
                .map(|(place, letters)| (offset + place, letters)),
        );
    }
    Some(joined)
}

/// The local ID that one reading gives its record at `place` (from 0): the counter starts at 1
/// and moves one per record (`FastaReader::generate_id`).
fn local_id(config: ReaderConfig, place: usize) -> String {
    format!("{}{}", config.id_prefix, place + 1)
}

/// One part: its place among the parts and its bytes `from..to`, and whether a part has failed.
struct Part<'a> {
    index: usize,
    from: u64,
    to: u64,
    failed: &'a std::sync::atomic::AtomicBool,
}

#[cfg(not(all(feature = "parallel", any(unix, target_os = "wasi"))))]
fn in_parts<T>(
    _input: Input<'_>,
    _start: u64,
    _threads: usize,
    _min_bytes: u64,
    _limits: Limits,
    _read: impl Fn(&mut At<'_>, &Part<'_>) -> Option<T>,
) -> Option<Vec<T>> {
    None
}

/// The results of `read` for the parts of `input` from byte `start` to its end, in order, read
/// by `threads` threads (each with its own `At`), or `None` where the reading is left as before:
/// fewer than two threads or `min_bytes` bytes, a part's start not found, or a part for which
/// `read` gives `None`.
#[cfg(all(feature = "parallel", any(unix, target_os = "wasi")))]
fn in_parts<T: Send>(
    input: Input<'_>,
    start: u64,
    threads: usize,
    min_bytes: u64,
    limits: Limits,
    read: impl Fn(&mut At<'_>, &Part<'_>) -> Option<T> + Sync,
) -> Option<Vec<T>> {
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    use std::sync::Mutex;

    let metadata = input.file.metadata().ok()?;
    if !metadata.is_file() {
        return None;
    }
    let end = metadata.len();
    let bytes = end.checked_sub(start)?;
    if threads < 2 || bytes == 0 || bytes < min_bytes {
        return None;
    }
    let parts = bytes
        .div_ceil(limits.max_span.max(1))
        .max(threads as u64 * limits.parts_per_thread.max(1));
    let span = bytes.div_ceil(parts);
    let parts = usize::try_from(bytes.div_ceil(span)).ok()?;
    let target = |part: usize| start + part as u64 * span;
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let results: Vec<Mutex<Option<T>>> = (0..parts).map(|_| Mutex::new(None)).collect();
    let work = || {
        let Some(mut at) = At::open(input) else {
            failed.store(true, Ordering::Relaxed);
            return;
        };
        while !failed.load(Ordering::Relaxed) {
            let part = next.fetch_add(1, Ordering::Relaxed);
            if part >= parts {
                return;
            }
            let from = if part == 0 {
                Some(start)
            } else {
                record_start(&mut at, start, target(part), end, limits.max_reach)
            };
            let to = if part + 1 == parts {
                Some(end)
            } else {
                record_start(&mut at, start, target(part + 1), end, limits.max_reach)
            };
            let result = match (from, to) {
                (Some(from), Some(to)) => read(
                    &mut at,
                    &Part {
                        index: part,
                        from,
                        to,
                        failed: &failed,
                    },
                ),
                _ => None,
            };
            match result {
                Some(result) => {
                    *results[part]
                        .lock()
                        .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(result);
                }
                None => failed.store(true, Ordering::Relaxed),
            }
        }
    };
    let helpers = threads.min(parts) - 1;
    std::thread::scope(|scope| {
        let spawned: Vec<_> = (0..helpers)
            .map(|_| std::thread::Builder::new().spawn_scoped(scope, work))
            .collect();
        work();
        for handle in spawned {
            match handle {
                Ok(handle) => {
                    if handle.join().is_err() {
                        failed.store(true, Ordering::Relaxed);
                    }
                }
                // A helper that could not start leaves its parts to the others.
                Err(_) => {}
            }
        }
    });
    if failed.load(Ordering::Relaxed) || next.load(Ordering::Relaxed) < parts {
        return None;
    }
    results
        .into_iter()
        .map(|part| {
            part.into_inner()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
        })
        .collect()
}

/// One thread's reads of the query file at offsets (where the parts read, the module's header).
#[cfg(unix)]
struct At<'a> {
    file: &'a File,
}

#[cfg(unix)]
impl<'a> At<'a> {
    fn open(input: Input<'a>) -> Option<At<'a>> {
        Some(At { file: input.file })
    }

    fn read_at(&mut self, buf: &mut [u8], at: u64) -> std::io::Result<usize> {
        std::os::unix::fs::FileExt::read_at(self.file, buf, at)
    }
}

/// One thread's reads of the query file at offsets (where the parts read, the module's header).
#[cfg(not(unix))]
struct At<'a> {
    file: File,
    /// The handle's position.
    position: u64,
    input: std::marker::PhantomData<Input<'a>>,
}

#[cfg(not(unix))]
impl<'a> At<'a> {
    /// `input`'s path opened again, where it is a regular file of `input.file`'s length and
    /// modification time.
    fn open(input: Input<'a>) -> Option<At<'a>> {
        let file = File::open(input.path).ok()?;
        let (first, again) = (input.file.metadata().ok()?, file.metadata().ok()?);
        let same = again.is_file()
            && again.len() == first.len()
            && again.modified().ok()? == first.modified().ok()?;
        same.then_some(At {
            file,
            position: 0,
            input: std::marker::PhantomData,
        })
    }

    fn read_at(&mut self, buf: &mut [u8], at: u64) -> std::io::Result<usize> {
        use std::io::{Read, Seek, SeekFrom};

        if self.position != at {
            self.position = self.file.seek(SeekFrom::Start(at))?;
        }
        let n = self.file.read(buf)?;
        self.position += n as u64;
        Ok(n)
    }
}

impl At<'_> {
    /// `buf.len()` bytes at `at`, or an error (an early end of the file among them).
    fn read_exact_at(&mut self, mut buf: &mut [u8], mut at: u64) -> std::io::Result<()> {
        while !buf.is_empty() {
            match self.read_at(buf, at) {
                Ok(0) => return Err(std::io::ErrorKind::UnexpectedEof.into()),
                Ok(n) => {
                    buf = &mut buf[n..];
                    at += n as u64;
                }
                Err(err) if err.kind() == std::io::ErrorKind::Interrupted => {}
                Err(err) => return Err(err),
            }
        }
        Ok(())
    }
}

/// The first record start at or after `target` (a `>` that begins a line after `start` and is not
/// followed by `?`), `end` if none comes before the end of the file, or `None` if none comes within
/// `max_reach` bytes after `target` before the end. Reads blocks of `BLOCK` bytes from `target`.
fn record_start(at: &mut At<'_>, start: u64, target: u64, end: u64, max_reach: u64) -> Option<u64> {
    const BLOCK: usize = 8192;
    if target >= end {
        return Some(end);
    }
    debug_assert!(target > start);
    // The record starts looked for are in `target..limit`.
    let limit = end.min(target.saturating_add(max_reach));
    let mut prev = [0u8];
    at.read_exact_at(&mut prev, target - 1).ok()?;
    let mut prev = prev[0];
    let mut block = vec![0u8; BLOCK];
    let mut position = target;
    while position < limit {
        let want = usize::try_from(end - position).map_or(BLOCK, |left| left.min(BLOCK));
        at.read_exact_at(&mut block[..want], position).ok()?;
        // A position is looked at with the byte after it in the block, or as the file's last byte.
        let at_end = position + want as u64 == end;
        let looked = if at_end { want } else { want - 1 };
        let looked = looked.min(usize::try_from(limit - position).unwrap_or(usize::MAX));
        for index in 0..looked {
            let byte = block[index];
            if byte == b'>' && prev == b'\n' && block[..want].get(index + 1) != Some(&b'?') {
                return Some(position + index as u64);
            }
            prev = byte;
        }
        position += looked as u64;
    }
    (limit == end).then_some(end)
}

/// What `keep` makes of each record of `part`, read by a new reader with the loop of
/// `x_query_stream::scan`, or `None` if its bytes hold a `\r` or a `>?`, cannot all be read, the
/// reading ends otherwise than at the end of the bytes, or another part failed.
fn read_part<T>(
    at: &mut At<'_>,
    part: &Part<'_>,
    config: ReaderConfig,
    keep: &mut dyn FnMut(FastaRecord) -> T,
) -> Option<Vec<T>> {
    let mut kept = Vec::new();
    let completed = read_part_with(at, part, config, false, &mut |record, _| {
        kept.push(keep(record));
        true
    })?;
    completed.then_some(kept)
}

/// The records of `part` as the reading before the search keeps them (with `pack` with their
/// packed letters), or `None` as for `read_part`, or if a record's ID or messages are not those of
/// one reading (the module's header).
fn read_part_records(
    at: &mut At<'_>,
    part: &Part<'_>,
    config: ReaderConfig,
    pack: bool,
) -> Option<PartRecords> {
    let mut kept = PartRecords::default();
    let completed = read_part_with(at, part, config, pack, &mut |record, packed| {
        if record.local_id != local_id(config, kept.records.len())
            || (part.index > 0 && !record.warnings.is_empty())
        {
            return false;
        }
        if let Some(letters) = packed {
            kept.packed.push((kept.records.len(), letters));
        }
        kept.records.push(record);
        true
    })?;
    completed.then_some(kept)
}

/// The loop of `x_query_stream::scan` over `part`'s bytes, each record given to `record` (with
/// `pack` with its packed letters, `FastaInputSource::x_pack_letters`) until it returns false;
/// `Some(true)` where every record was taken and the reading ended at the end of the bytes without
/// a `\r` or a `>?` in them, `Some(false)` where `record` refused one, `None` otherwise.
fn read_part_with(
    at: &mut At<'_>,
    part: &Part<'_>,
    config: ReaderConfig,
    pack: bool,
    record: &mut dyn FnMut(FastaRecord, Option<PackedLetters>) -> bool,
) -> Option<bool> {
    use crate::blastinput::fasta_reader::FastaInputSource;
    use std::sync::atomic::Ordering;

    if part.from >= part.to {
        return Some(true);
    }
    let mut refused = false;
    let mut source = FastaInputSource::x_from_reader(
        PartReader {
            at,
            position: part.from,
            end: part.to,
            last: b'\n',
            refused: &mut refused,
        },
        config,
    );
    if pack {
        source.x_pack_letters();
    }
    loop {
        if source.end() {
            break;
        }
        if part.failed.load(Ordering::Relaxed) {
            return None;
        }
        match source.next_sequence(&mut |_| Ok(())) {
            Ok(read) => {
                let packed = if pack { source.x_take_packed() } else { None };
                if !record(read, packed) {
                    return Some(false);
                }
            }
            Err(_) => return None,
        }
    }
    drop(source);
    (!refused).then_some(true)
}

/// Bytes `position..end` of a file, read at their offsets (`At`), that notes (`refused`) a `\r`
/// or a `>?` in what it returns, and an I/O error or an end of the file before `end` (which the
/// reader would take for the end of its input).
struct PartReader<'a, 'b> {
    at: &'a mut At<'b>,
    position: u64,
    end: u64,
    /// The last byte returned (a part starts at the start of a line).
    last: u8,
    refused: &'a mut bool,
}

impl std::io::Read for PartReader<'_, '_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let left = usize::try_from(self.end - self.position).unwrap_or(usize::MAX);
        let want = buf.len().min(left);
        if want == 0 {
            return Ok(0);
        }
        let n = match self.at.read_at(&mut buf[..want], self.position) {
            Ok(0) | Err(_) => {
                *self.refused = true;
                return Ok(0);
            }
            Ok(n) => n,
        };
        let bytes = &buf[..n];
        if bytes.contains(&b'\r')
            || (bytes.contains(&b'?')
                && (bytes[0] == b'?' && self.last == b'>'
                    || bytes.windows(2).any(|pair| pair == b">?")))
        {
            *self.refused = true;
        }
        self.last = bytes[n - 1];
        self.position += n as u64;
        Ok(n)
    }
}

#[cfg(all(test, unix, feature = "parallel"))]
mod tests {
    use super::*;
    use crate::blastinput::fasta_reader::{FastaInputSource, ParseErrorCode, ReadError};
    use std::io::Write;

    /// The first reading as `x_query_stream::scan` makes it (without the archive): the lengths,
    /// and whether it ended at the end of the input.
    fn serial(bytes: &[u8], config: ReaderConfig) -> (Vec<usize>, bool) {
        let mut source = FastaInputSource::from_bytes(bytes, config);
        let mut lengths = Vec::new();
        loop {
            if source.end() {
                return (lengths, true);
            }
            match source.next_sequence(&mut |_| Ok(())) {
                Ok(record) => lengths.push(record.sequence.len()),
                Err(ReadError::Parse {
                    code: ParseErrorCode::Eof,
                    ..
                }) => return (lengths, false),
                Err(_) => return (lengths, false),
            }
        }
    }

    fn file_of(bytes: &[u8], tag: &str) -> (std::path::PathBuf, File) {
        let path = std::env::temp_dir().join(format!(
            "losat_qscanpar_{}_{tag}_{:?}.fa",
            std::process::id(),
            std::thread::current().id()
        ));
        let mut file = File::create(&path).expect("create the test file");
        file.write_all(bytes).expect("write the test file");
        drop(file);
        let file = File::open(&path).expect("open the test file");
        (path, file)
    }

    fn at_of<'a>(file: &'a File, path: &'a std::path::Path) -> At<'a> {
        At::open(Input { file, path }).expect("the test file")
    }

    /// The parallel reading of `bytes` from byte `start`, without a minimum size.
    fn parallel(
        bytes: &[u8],
        start: u64,
        config: ReaderConfig,
        threads: usize,
        limits: Limits,
    ) -> Option<Vec<usize>> {
        let (path, file) = file_of(bytes, "p");
        let input = Input {
            file: &file,
            path: &path,
        };
        let lengths = lengths_with(input, start, config, threads, 0, limits);
        let _ = std::fs::remove_file(path);
        lengths
    }

    /// A small deterministic generator (xorshift).
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    /// A FASTA text with the reader's cases: ids and none, titles, lower case, IUPAC letters, `N`
    /// runs, digits and spaces in data lines, blank and comment lines, deflines without letters,
    /// data before the first defline, lines of any length and no final newline; in one text of four
    /// `>?` gap lines and `>?_` deflines, in one of eight a line of bad letters (a reader error).
    fn random_fasta(rng: &mut Rng, records: usize) -> Vec<u8> {
        let mut text = Vec::new();
        let gaps = rng.below(4) == 0;
        let bad = rng.below(8) == 0;
        if rng.below(4) == 0 {
            text.extend_from_slice(b"ACGTNACGT\n");
        }
        for index in 0..records {
            match rng.below(6) {
                0 => text.extend_from_slice(b">\n"),
                1 => text.extend_from_slice(format!(">lcl|q{index} title {index}\n").as_bytes()),
                2 => text.extend_from_slice(format!(">q{index}\n").as_bytes()),
                3 if gaps => text.extend_from_slice(b">?_x a title\n"),
                _ => text.extend_from_slice(format!(">gnl|db|{index} some words\n").as_bytes()),
            }
            let lines = rng.below(8);
            for _ in 0..lines {
                match rng.below(14) {
                    0 => text.push(b'\n'),
                    1 => text.extend_from_slice(b"; a comment\n"),
                    2 if gaps => text.extend_from_slice(b">?10\n"),
                    3 => text.extend_from_slice(b"  12 acgt acgt 13\n"),
                    4 if bad && rng.below(4) == 0 => text.extend_from_slice(b"%%%%%%%%%%%%\n"),
                    _ => {
                        let width = 1 + rng.below(90) as usize;
                        for _ in 0..width {
                            let letters = b"ACGTacgtNNnRYKMSWBDHV";
                            text.push(letters[rng.below(letters.len() as u64) as usize]);
                        }
                        text.push(b'\n');
                    }
                }
            }
        }
        if !text.is_empty() && rng.below(3) == 0 {
            text.pop();
        }
        text
    }

    const SMALL: Limits = Limits {
        max_span: 97,
        max_reach: 400,
        parts_per_thread: 3,
    };

    #[test]
    fn the_parts_give_the_first_readings_lengths() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let (mut taken, mut refused, mut failed) = (0, 0, 0);
        for round in 0..300 {
            let loaders = round % 2 == 0;
            let config = ReaderConfig::query("BLASTN", false, loaders);
            let text = random_fasta(&mut rng, 1 + round % 40);
            let (expected, clean) = serial(&text, config);
            for threads in [2, 3, 8] {
                for limits in [
                    SMALL,
                    Limits {
                        max_span: 7,
                        max_reach: 2000,
                        parts_per_thread: 1,
                    },
                    Limits::DEFAULT,
                ] {
                    match parallel(&text, 0, config, threads, limits) {
                        Some(lengths) => {
                            assert!(clean, "round {round}: parts taken where the reading fails");
                            assert!(!has_gap_mark(&text), "round {round}: parts taken with `>?`");
                            assert_eq!(lengths, expected, "round {round} threads {threads}");
                            taken += 1;
                        }
                        None => refused += 1,
                    }
                    if !clean {
                        failed += 1;
                    }
                }
            }
        }
        assert!(taken > 1000, "the parts were taken {taken} times");
        assert!(
            refused > 300 && failed > 30,
            "refused {refused}, failed readings {failed}"
        );
    }

    fn has_gap_mark(text: &[u8]) -> bool {
        text.windows(2).any(|pair| pair == b">?")
    }

    /// The record start of `record_start` found by looking at every byte.
    fn record_start_by_bytes(text: &[u8], target: usize, reach: usize) -> Option<usize> {
        if target >= text.len() {
            return Some(text.len());
        }
        let limit = text.len().min(target + reach);
        (target..limit)
            .find(|&at| {
                text[at] == b'>' && text[at - 1] == b'\n' && text.get(at + 1) != Some(&b'?')
            })
            .or((limit == text.len()).then_some(text.len()))
    }

    #[test]
    fn record_starts_across_blocks_are_those_of_every_byte() {
        let mut rng = Rng(0x2545_f491_4f6c_dd1d);
        for round in 0..6 {
            let mut text = random_fasta(&mut rng, 120 + 40 * round);
            // Record starts at and around the block edges of the search (8,192 bytes from a target).
            for edge in [8_190usize, 8_191, 8_192, 8_193, 16_383, 16_384] {
                if edge + 2 < text.len() {
                    text[edge - 1] = b'\n';
                    text[edge] = b'>';
                }
            }
            let (path, file) = file_of(&text, "rsb");
            let mut at = at_of(&file, &path);
            let end = text.len() as u64;
            for target in (1..text.len())
                .step_by(37)
                .chain([1, 2, 8_190, 8_191, 8_192, 8_193])
            {
                for reach in [1usize, 50, 8_191, 8_192, 8_193, 20_000, 1 << 20] {
                    assert_eq!(
                        record_start(&mut at, 0, target as u64, end, reach as u64),
                        record_start_by_bytes(&text, target, reach).map(|at| at as u64),
                        "round {round} target {target} reach {reach}"
                    );
                }
            }
            let _ = std::fs::remove_file(path);
        }
    }

    #[test]
    fn a_part_refuses_gap_marks_across_reads_and_bytes_it_cannot_read() {
        let config = ReaderConfig::query("BLASTN", false, true);
        let failed = std::sync::atomic::AtomicBool::new(false);
        // The reader's input window is 8,191 bytes: a `>` as its last byte, `?` the next read's first.
        let mut text = b">a\n".to_vec();
        while text.len() < 8_189 {
            text.extend_from_slice(b"ACGTACGTAC\n");
        }
        text.truncate(8_189);
        text.push(b'\n');
        let mut gap = text.clone();
        gap.extend_from_slice(b">?5\nACGT\n");
        assert_eq!(gap[8_190..8_192], *b">?");
        let mut plain = text.clone();
        plain.extend_from_slice(b">b\nACGT\n");
        for (bytes, taken) in [(&gap, false), (&plain, true)] {
            let (path, file) = file_of(bytes, "gap");
            let mut at = at_of(&file, &path);
            let end = bytes.len() as u64;
            let lengths = |at: &mut At<'_>, to: u64| {
                let part = Part {
                    index: 0,
                    from: 0,
                    to,
                    failed: &failed,
                };
                read_part(at, &part, config, &mut |record| record.sequence.len())
            };
            assert_eq!(lengths(&mut at, end).is_some(), taken);
            // A part past the end of the file: the reader would take the early end for its own.
            assert_eq!(lengths(&mut at, end + 100), None);
            let _ = std::fs::remove_file(path);
        }
    }

    #[test]
    fn a_reading_that_fails_or_cannot_be_cut_is_left_to_the_first_reading() {
        let config = ReaderConfig::query("BLASTN", false, true);
        let good = b">a\nACGT\n>b\nACGTACGT\n>c\nAC\n>d\nACGTA\n".to_vec();
        assert_eq!(parallel(&good, 0, config, 4, SMALL), Some(vec![4, 8, 2, 5]));
        // A `\r` anywhere.
        let mut cr = good.clone();
        cr.splice(20..20, *b"\r");
        assert_eq!(parallel(&cr, 0, config, 4, SMALL), None);
        // A reader error in a late part (a line of mostly bad letters).
        let mut bad = good.clone();
        bad.extend_from_slice(b">e\n%%%%%%%%%%%%\n");
        assert!(!serial(&bad, config).1);
        assert_eq!(parallel(&bad, 0, config, 4, SMALL), None);
        // A record longer than the reach of a target.
        let mut long = b">a\n".to_vec();
        long.extend(std::iter::repeat_n(b'A', 5000));
        long.extend_from_slice(b"\n>b\nAC\n");
        let limits = Limits {
            max_span: 100,
            max_reach: 200,
            parts_per_thread: 1,
        };
        assert_eq!(parallel(&long, 0, config, 4, limits), None);
        // A Seq-id line first (the data loaders' check of the first line), and comment lines only
        // (an `eEOF`): the parts are taken only where the first reading ends at the end.
        for text in [
            b"NC_000913.3\n>a\nACGT\n>b\nAC\n".as_slice(),
            b"; one\n; two\n\n",
        ] {
            let (expected, clean) = serial(text, config);
            let parts = parallel(text, 0, config, 4, SMALL);
            assert_eq!(
                parts.is_some(),
                clean,
                "{:?}",
                String::from_utf8_lossy(text)
            );
            if let Some(lengths) = parts {
                assert_eq!(lengths, expected);
            }
        }
        assert!(!serial(b"; one\n; two\n\n", config).1);
        // One thread, or fewer bytes than the minimum.
        assert_eq!(parallel(&good, 0, config, 1, SMALL), None);
        let (path, file) = file_of(&good, "min");
        let input = Input {
            file: &file,
            path: &path,
        };
        assert_eq!(lengths_with(input, 0, config, 4, 1 << 20, SMALL), None);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn the_parts_start_where_the_reader_starts() {
        let config = ReaderConfig::query("BLASTN", false, true);
        let mut text = b"   \n\n".to_vec();
        let body = b">a\nACGT\n>b\nACGTACGT\n>c\nAC\n>d\nACGTA\n";
        text.extend_from_slice(body);
        let start = 5;
        assert_eq!(
            parallel(&text, start, config, 3, SMALL),
            Some(serial(&text[start as usize..], config).0)
        );
    }

    #[test]
    fn record_starts_skip_gap_lines_and_lines_that_do_not_begin_with_the_mark() {
        let text = b">a\nAC\n>?5\n x>y\n>b\nAC\n".to_vec();
        let (path, file) = file_of(&text, "rs");
        let mut at = at_of(&file, &path);
        let end = text.len() as u64;
        assert_eq!(record_start(&mut at, 0, 1, end, 100), Some(15));
        assert_eq!(record_start(&mut at, 0, 15, end, 100), Some(15));
        assert_eq!(record_start(&mut at, 0, 16, end, 100), Some(end));
        assert_eq!(record_start(&mut at, 0, 1, end, 5), None);
        assert_eq!(record_start(&mut at, 0, 1, end, 15), Some(15));
        let _ = std::fs::remove_file(path);
    }

    /// The reading before the search as `search_cli` makes it without the parts: the records (with
    /// `pack` their packed letters), and whether it ended at the end of the input.
    fn serial_records(bytes: &[u8], config: ReaderConfig, pack: bool) -> (PartRecords, bool) {
        use crate::blastinput::fasta_reader::{
            read_queries, x_read_queries_packed, QueryEnd, QueryRecords,
        };
        let mut source = FastaInputSource::from_bytes(bytes, config);
        let (QueryRecords { records, end }, packed) = if pack {
            x_read_queries_packed(&mut source)
        } else {
            (read_queries(&mut source), Vec::new())
        };
        (
            PartRecords { records, packed },
            matches!(end, QueryEnd::Input),
        )
    }

    /// `records` of `bytes` from byte `start`, without a minimum size.
    fn parallel_records(
        bytes: &[u8],
        start: u64,
        config: ReaderConfig,
        threads: usize,
        pack: bool,
        limits: Limits,
    ) -> Option<PartRecords> {
        let (path, file) = file_of(bytes, "r");
        let input = Input {
            file: &file,
            path: &path,
        };
        let records = records_with(input, start, config, threads, pack, 0, limits);
        let _ = std::fs::remove_file(path);
        records
    }

    /// A FASTA text of records with titles and letters only (upper and lower case, `N` runs, lines
    /// of any length), of which the reader writes no message.
    fn plain_fasta(rng: &mut Rng, records: usize) -> Vec<u8> {
        let mut text = Vec::new();
        for index in 0..records {
            text.extend_from_slice(format!(">q{index} a title {index}\n").as_bytes());
            for _ in 0..1 + rng.below(6) {
                let width = 1 + rng.below(90) as usize;
                for _ in 0..width {
                    let letters = b"ACGTACGTacgtNn";
                    text.push(letters[rng.below(letters.len() as u64) as usize]);
                }
                text.push(b'\n');
            }
        }
        text
    }

    #[test]
    fn the_parts_give_the_readings_records() {
        let mut rng = Rng(0x6a09_e667_f3bc_c908);
        let (mut taken, mut refused) = (0, 0);
        for round in 0..200 {
            let config = ReaderConfig::query("BLASTN", false, round % 4 < 2);
            // Every other text of records the reader writes no message for, as most query files.
            let mut text = if round % 2 == 0 {
                random_fasta(&mut rng, 1 + round % 40)
            } else {
                plain_fasta(&mut rng, 1 + round % 40)
            };
            // In one text of five a hyphen in a data line (a message of the reader).
            if round % 5 == 0 {
                if let Some(at) = text.windows(3).position(|bytes| bytes == b"\nAC") {
                    text.insert(at + 2, b'-');
                }
            }
            for pack in [false, true] {
                let (expected, clean) = serial_records(&text, config, pack);
                for threads in [2, 8] {
                    for limits in [
                        SMALL,
                        Limits {
                            max_span: 7,
                            max_reach: 2000,
                            parts_per_thread: 1,
                        },
                    ] {
                        match parallel_records(&text, 0, config, threads, pack, limits) {
                            Some(parts) => {
                                assert!(
                                    clean,
                                    "round {round}: parts taken where the reading fails"
                                );
                                assert_eq!(parts.records, expected.records, "round {round}");
                                assert_eq!(parts.packed, expected.packed, "round {round}");
                                taken += 1;
                            }
                            None => refused += 1,
                        }
                    }
                }
            }
        }
        assert!(
            taken > 300 && refused > 100,
            "taken {taken}, refused {refused}"
        );
    }

    #[test]
    fn records_with_messages_after_the_first_part_are_left_to_the_reading() {
        let config = ReaderConfig::query("BLASTN", false, true);
        let late = b">a\nACGT\n>b\nACGTACGT\n>c\nAC\n>d\nAC-GTA\n";
        let (expected, clean) = serial_records(late, config, false);
        assert!(clean && !expected.records[3].warnings.is_empty());
        assert!(parallel_records(late, 0, config, 4, false, SMALL).is_none());
        // The same record first: its messages count lines as one reading does.
        let first = b">d\nAC-GTA\n>a\nACGT\n>b\nACGTACGT\n>c\nAC\n";
        let (expected, clean) = serial_records(first, config, false);
        assert!(clean && !expected.records[0].warnings.is_empty());
        let parts = parallel_records(first, 0, config, 4, false, SMALL).expect("the parts");
        assert_eq!(parts.records, expected.records);
        let ids: Vec<&str> = parts
            .records
            .iter()
            .map(|record| record.local_id.as_str())
            .collect();
        assert_eq!(ids, ["Query_1", "Query_2", "Query_3", "Query_4"]);
    }
}
