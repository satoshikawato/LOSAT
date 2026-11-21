//! The first reading of LOSAT_X_QSTREAM (`x_query_stream::scan`, which keeps each query record's
//! length and how the reading ends) cut into parts at record starts and read on `-num_threads`
//! threads (the calling thread and scoped helpers, before the search pool exists), each part by a
//! reader of its own. It is the form of that reading wherever LOSAT_X_QSTREAM acts on a file of
//! `MIN_BYTES` or more: it had its own switch, LOSAT_X_QSCANPAR, until 2026-11-01, when the switch
//! was folded into LOSAT_X_QSTREAM (line SF S3, family bn-query-read). The minimum stays: forced
//! below it, blastn files of 2.1 and 4.0 MB at two threads on the workstation's idle socket ran
//! 0.07 to 0.32 s longer in 11 of 32 runs (against 1 of 32 with the one-thread reading); the cause
//! is not known.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
//! NCBI reads its query batches one at a time, inside the loop. LOSAT needs every record's length
//! before its first batch (the batch planner, LOSAT_X_RPTSTREAM's checks),
//! so where LOSAT_X_QSTREAM acts it reads the whole query file once for the lengths alone, on one
//! thread, while every other thread of the search waits: for dc-megablast of 31,168 D. melanogaster
//! cDNAs (96 MB) against 3 Mb at 16 threads on the workstation, 0.12 s of a 1.55 s search. No NCBI
//! counterpart: it is speed only.
//!
//! What it does. The bytes from the reader's start to the end of the file are cut at `parts`
//! targets (at most `MAX_SPAN` bytes apart, at least four per thread). Each part runs from the
//! first record start at or after its target (a `>` at the start of a line, not followed by `?`)
//! to the next part's. The threads take the parts in turn; each part is read by a new reader over
//! its bytes (`FastaInputSource::x_from_reader`, the same configuration) with the loop of `scan`:
//! `end()`, then `next_sequence`, the length of each record kept and the record dropped. The
//! lengths of the parts, in order, are the result.
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
//! When it does not act (the first reading runs on one thread, `scan` on the search's reader, which
//! no part has moved: the parts read with `read_at`, which leaves the file position alone): one
//! thread; a target without `unix` file reads at an offset or without the `parallel` feature
//! (WebAssembly, Windows); the BLAST archive (the first reading notes the titles it cannot write,
//! `x_query_stream::scan`); a file of fewer than `MIN_BYTES` bytes after the reader's start
//! (`LOSAT_X_QSCANPARMIN`, tuning, sets it); and, found while the parts are read, a `\r` anywhere in
//! the file, a `>?` anywhere in it (a gap line, whose length becomes that many `N`s, or a `>?_`
//! defline), a target with no record start within `MAX_REACH` bytes after it (a record that long
//! would take that much memory in every thread that reads one), a part whose bytes cannot all be
//! read (an I/O error or a file that ends early: a reader takes either for the end of its input),
//! or a part whose reading ends with anything but the end of its bytes (an `eEOF`, a reader error or
//! LOSAT's rejection of a Seq-id line: the first reading then finds it, with its messages, as
//! before). What is read twice then is the parts' bytes up to the failure.
//!
//! Memory: a part holds at most `MAX_SPAN + MAX_REACH` bytes of records, and without gap lines a
//! record has no more letters than bytes, so a thread holds one record of at most that many letters
//! and the reader's buffers at a time. No floating-point operation is involved.
//!
//! With `LOSAT_X_STATS`, `scan_par=1` in LOSAT_X_QSTREAM's line says that the parts' lengths were
//! taken (`scan_par=0`: the one-thread reading).
//!
//! The second reading (since 2026-11-03). Where the parts' lengths were taken, LOSAT_X_QSTREAM's
//! second reading reads the same parts again (`Span`: each part's bytes, first record, records, the
//! lines before it and its lines, the line number of a reader that read it to its end), when the
//! batch loop needs their records, on the search's threads (`x_query_stream::State::read_parts`).
//! Each part is read by a new reader that starts with the line number and the local-ID counter the
//! one reader has at the part's first byte (`read_again`, `FastaInputSource::x_from_reader_at`):
//! the two values above are all that a record's reading takes from the records before it, so the
//! part's records are the one reader's, their local IDs and the line numbers of their messages
//! included (so, unlike the parts of the store, a record with a message in a later part needs no
//! one reader). `read_again` checks that the file did not change in a way that would make that
//! untrue: the part starts at a record start, its bytes end where the first reading's did and hold
//! no `\r` or `>?`, its records have the first reading's number and lengths, and its lines are the
//! first reading's; the last part is read to the file's end, up to its records, as the one reader
//! reads it.
//!
//! The records of the store (EXPERIMENT LOSAT_X_BNQPARTS, `store_parts`). Where LOSAT_X_BNQPACK's
//! store holds the records read before the search, the same parts of the file (the same cut, the
//! same `read_part` conditions) are read on the same threads into the parts of the store
//! (`x_query_drop::Part`): each part by a new reader that packs the letters as the search's reader
//! does, every record kept, in order. A part's records are the first reading's records of those
//! bytes (the argument above), but for what a new reader numbers from its own start: the local ID
//! (its counter starts at 1 in each part; the store gives each record the ID of its position among
//! all records, so a record's ID is checked to be its part reader's counter's) and the line of a
//! reader message, so a record of a part after the first that has a message leaves the file to the
//! search's reader. The parts after the first need the parse of their records to have no state
//! across the file: `-parse_deflines` (its records with ids take no number of the counter, and its
//! reader adds every record to the search's scope of ids) and the archive are not read in parts.

#![cfg_attr(not(all(unix, feature = "parallel")), allow(dead_code))]

use std::fs::File;
use std::ops::Range;
use std::sync::OnceLock;

use super::x_query_drop::Part;
#[cfg(all(unix, feature = "parallel"))]
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::fasta_reader::ReaderConfig;

/// The bytes after the reader's start below which the first reading stays on one thread: about
/// 4 ms of it on the workstation.
const MIN_BYTES: u64 = 4 << 20;
/// The most bytes between two targets.
const MAX_SPAN: u64 = 256 << 10;
/// The most bytes after a target searched for a record start.
const MAX_REACH: u64 = 256 << 10;
/// The parts per thread, at least.
const PARTS_PER_THREAD: u64 = 4;

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

/// A part of the first reading that holds records, for the second reading (`read_again`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Span {
    /// The part's bytes in the file.
    pub(super) bytes: Range<u64>,
    /// Its first record (the input record number) and its records.
    pub(super) first: usize,
    pub(super) records: usize,
    /// The lines before it and its lines (the line number of a reader that read it to its end).
    pub(super) line: u64,
    pub(super) lines: u64,
}

/// The first reading in parts: every record's length, and the parts that hold records.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct Scan {
    pub(super) lengths: Vec<usize>,
    pub(super) spans: Vec<Span>,
}

/// The lengths of the query records from byte `start` of `file` to its end, read in parts by
/// `threads` threads, with the parts, or `None` where the first reading runs as before (the
/// module's header). `config` is the first reading's configuration, without the BLAST archive.
pub(super) fn lengths(
    file: &File,
    start: u64,
    config: ReaderConfig,
    threads: usize,
) -> Option<Scan> {
    debug_assert!(!config.keep_defline);
    lengths_with(file, start, config, threads, min_bytes(), Limits::DEFAULT)
}

/// `lengths` without a minimum size and with small parts, for the tests of the second reading
/// (`x_query_stream.rs`).
#[cfg(all(test, unix, feature = "parallel"))]
pub(super) fn small_lengths(
    file: &File,
    start: u64,
    config: ReaderConfig,
    threads: usize,
) -> Option<Scan> {
    let limits = Limits {
        max_span: 4 << 10,
        max_reach: 2 << 20,
        parts_per_thread: 4,
    };
    lengths_with(file, start, config, threads, 0, limits)
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

#[cfg(not(all(unix, feature = "parallel")))]
fn lengths_with(
    _file: &File,
    _start: u64,
    _config: ReaderConfig,
    _threads: usize,
    _min_bytes: u64,
    _limits: Limits,
) -> Option<Scan> {
    None
}

#[cfg(all(unix, feature = "parallel"))]
fn lengths_with(
    file: &File,
    start: u64,
    config: ReaderConfig,
    threads: usize,
    min_bytes: u64,
    limits: Limits,
) -> Option<Scan> {
    let parts = parts_with(
        file,
        start,
        threads,
        min_bytes,
        limits,
        |from, to, _, failed| {
            let (lengths, lines) = read_part(file, from..to, config, failed)?;
            Some((from..to, lengths, lines))
        },
    )?;
    let mut scan = Scan {
        lengths: Vec::new(),
        spans: Vec::new(),
    };
    let mut line = 0;
    for (bytes, lengths, lines) in parts {
        // An empty part (two targets found the same record start) holds no records.
        if !lengths.is_empty() {
            scan.spans.push(Span {
                bytes,
                first: scan.lengths.len(),
                records: lengths.len(),
                line,
                lines,
            });
        }
        line += lines;
        scan.lengths.extend(lengths);
    }
    Some(scan)
}

/// The results of `read` on the parts of `file` from byte `start` to its end (the module's header),
/// in order, read by `threads` threads, or `None` where a part could not be cut or read (`read`
/// gives `None`; it is called with the part's bytes `from..to`, its index and a flag set when
/// another part failed).
#[cfg(all(unix, feature = "parallel"))]
fn parts_with<T: Send>(
    file: &File,
    start: u64,
    threads: usize,
    min_bytes: u64,
    limits: Limits,
    read: impl Fn(u64, u64, usize, &std::sync::atomic::AtomicBool) -> Option<T> + Sync,
) -> Option<Vec<T>> {
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    use std::sync::Mutex;

    let end = file.metadata().ok()?.len();
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
        while !failed.load(Ordering::Relaxed) {
            let part = next.fetch_add(1, Ordering::Relaxed);
            if part >= parts {
                return;
            }
            let from = if part == 0 {
                Some(start)
            } else {
                record_start(file, start, target(part), end, limits.max_reach)
            };
            let to = if part + 1 == parts {
                Some(end)
            } else {
                record_start(file, start, target(part + 1), end, limits.max_reach)
            };
            let got = match (from, to) {
                (Some(from), Some(to)) => read(from, to, part, &failed),
                _ => None,
            };
            match got {
                Some(got) => {
                    *results[part]
                        .lock()
                        .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(got);
                }
                None => failed.store(true, Ordering::Relaxed),
            }
        }
    };
    let helpers = threads.min(parts) - 1;
    on_helpers(helpers, &failed, work);
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

/// Runs `work` on the calling thread and on `helpers` scoped threads, and waits for them; a helper
/// that panics sets `failed`, one that cannot start leaves its work to the others.
#[cfg(all(unix, feature = "parallel"))]
fn on_helpers(helpers: usize, failed: &std::sync::atomic::AtomicBool, work: impl Fn() + Sync) {
    let work = &work;
    std::thread::scope(|scope| {
        let spawned: Vec<_> = (0..helpers)
            .map(|_| std::thread::Builder::new().spawn_scoped(scope, work))
            .collect();
        work();
        // A helper that could not start leaves its parts to the others.
        for handle in spawned.into_iter().flatten() {
            if handle.join().is_err() {
                failed.store(true, std::sync::atomic::Ordering::Relaxed);
            }
        }
    });
}

/// EXPERIMENT (LOSAT_X_BNQPARTS): the query records from byte `start` of `file` to its end, read in
/// parts by `threads` threads into the parts of the store (`x_query_drop::Part`, their letters
/// packed as the reader packs them), or `None` where the search's reader reads them into one part
/// (the module's header, "The records of the store"). `config` is the reader's configuration,
/// without the BLAST archive or `-parse_deflines`.
pub(super) fn store_parts(
    file: &File,
    start: u64,
    config: ReaderConfig,
    threads: usize,
) -> Option<Vec<Part>> {
    debug_assert!(!config.keep_defline && !config.parse_deflines);
    store_parts_with(file, start, config, threads, min_bytes(), Limits::DEFAULT)
}

#[cfg(not(all(unix, feature = "parallel")))]
fn store_parts_with(
    _file: &File,
    _start: u64,
    _config: ReaderConfig,
    _threads: usize,
    _min_bytes: u64,
    _limits: Limits,
) -> Option<Vec<Part>> {
    None
}

#[cfg(all(unix, feature = "parallel"))]
fn store_parts_with(
    file: &File,
    start: u64,
    config: ReaderConfig,
    threads: usize,
    min_bytes: u64,
    limits: Limits,
) -> Option<Vec<Part>> {
    parts_with(
        file,
        start,
        threads,
        min_bytes,
        limits,
        |from, to, part, failed| read_part_store(file, from, to, config, part, failed),
    )
}

/// The records of bytes `from..to` of `file` (part `part`), read by a new reader as
/// `x_read_queries_packed` reads them into a part of the store, or `None` where `read_part` gives
/// `None`, or a record of a part after the first has a reader message (its text may name a line,
/// counted from the part's start), or a record's local ID is not the reader's counter's.
#[cfg(all(unix, feature = "parallel"))]
fn read_part_store(
    file: &File,
    from: u64,
    to: u64,
    config: ReaderConfig,
    part: usize,
    failed: &std::sync::atomic::AtomicBool,
) -> Option<Part> {
    // The store of a part takes about a third of its bytes (titles, and a quarter of a byte per
    // letter); it grows where it needs more.
    let span = usize::try_from(to.saturating_sub(from)).ok()?;
    let mut store = Part::with_capacity(span / 2, span / 128);
    read_records(
        file,
        from..to,
        config,
        (0, 1),
        true,
        failed,
        |record, packed| {
            if part > 0 && !record.warnings.is_empty() {
                return After::Refuse;
            }
            let n = store.len() + 1;
            if store.push(record, packed, config.id_prefix, n) {
                After::Next
            } else {
                After::Refuse
            }
        },
    )
    .ok()?;
    store.shrink();
    Some(store)
}

/// The first record start at or after `target` (a `>` that begins a line after `start` and is not
/// followed by `?`), `end` if none comes before the end of the file, or `None` if none comes within
/// `max_reach` bytes after `target` before the end. Reads blocks of `BLOCK` bytes from `target`.
#[cfg(all(unix, feature = "parallel"))]
fn record_start(file: &File, start: u64, target: u64, end: u64, max_reach: u64) -> Option<u64> {
    use std::os::unix::fs::FileExt;

    const BLOCK: usize = 8192;
    if target >= end {
        return Some(end);
    }
    debug_assert!(target > start);
    // The record starts looked for are in `target..limit`.
    let limit = end.min(target.saturating_add(max_reach));
    let mut prev = [0u8];
    file.read_exact_at(&mut prev, target - 1).ok()?;
    let mut prev = prev[0];
    let mut block = vec![0u8; BLOCK];
    let mut at = target;
    while at < limit {
        let want = usize::try_from(end - at).map_or(BLOCK, |left| left.min(BLOCK));
        file.read_exact_at(&mut block[..want], at).ok()?;
        // A position is looked at with the byte after it in the block, or as the file's last byte.
        let at_end = at + want as u64 == end;
        let looked = if at_end { want } else { want - 1 };
        let looked = looked.min(usize::try_from(limit - at).unwrap_or(usize::MAX));
        for index in 0..looked {
            let byte = block[index];
            if byte == b'>' && prev == b'\n' && block[..want].get(index + 1) != Some(&b'?') {
                return Some(at + index as u64);
            }
            prev = byte;
        }
        at += looked as u64;
    }
    (limit == end).then_some(end)
}

/// What a part's reading does after a record (`read_records`).
#[cfg(all(unix, feature = "parallel"))]
enum After {
    /// Reads the next record, if its bytes hold one.
    Next,
    /// Stops: the records wanted are read.
    Stop,
    /// Gives the part up.
    Refuse,
}

/// Why a part's reading failed (`read_records`).
#[cfg(all(unix, feature = "parallel"))]
enum PartError {
    /// Its bytes hold a `\r` or a `>?`, or could not all be read (`PartReader`).
    Bytes,
    /// The reader failed: an `eEOF`, a reader error or LOSAT's rejection of a Seq-id line.
    Reader(crate::blastinput::fasta_reader::ReadError),
    /// A record was refused (`After::Refuse`), or another part failed (`failed`).
    Refused,
}

/// Reads the records of bytes `bytes` of `file` with a new reader that starts at `start` (the lines
/// before it and the local ID of its first record, `x_from_reader_at`; `(0, 1)` is a new reader's)
/// with the loop of `x_query_stream::scan` (`end()`, then `next_sequence`), the letters packed as
/// the search's reader packs them where `pack` is set, and gives each record (and its packed
/// letters) to `each`. Returns the lines read (those of the part), or why the reading failed; it
/// stops when another part failed (`failed`).
#[cfg(all(unix, feature = "parallel"))]
fn read_records(
    file: &File,
    bytes: Range<u64>,
    config: ReaderConfig,
    start: (u64, i64),
    pack: bool,
    failed: &std::sync::atomic::AtomicBool,
    mut each: impl FnMut(FastaRecord, Option<crate::sequence::packed_letters::PackedLetters>) -> After,
) -> Result<u64, PartError> {
    use crate::blastinput::fasta_reader::FastaInputSource;
    use std::sync::atomic::Ordering;

    if bytes.start >= bytes.end {
        return Ok(0);
    }
    let mut refused = false;
    let mut source = FastaInputSource::x_from_reader_at(
        PartReader {
            file,
            at: bytes.start,
            end: bytes.end,
            last: b'\n',
            refused: &mut refused,
        },
        config,
        start.0,
        start.1,
    );
    if pack {
        source.x_pack_letters();
    }
    loop {
        if source.end() {
            break;
        }
        if failed.load(Ordering::Relaxed) {
            return Err(PartError::Refused);
        }
        match source.next_sequence(&mut |_| Ok(())) {
            Ok(record) => {
                let packed = if pack { source.x_take_packed() } else { None };
                match each(record, packed) {
                    After::Next => {}
                    After::Stop => break,
                    After::Refuse => return Err(PartError::Refused),
                }
            }
            Err(error) => return Err(PartError::Reader(error)),
        }
    }
    // The reader counts from the lines before the part.
    let lines = source.x_line_number() - start.0;
    drop(source);
    if refused {
        return Err(PartError::Bytes);
    }
    Ok(lines)
}

/// The lengths of the records of bytes `bytes` of `file` and the lines read, as the first reading
/// reads a part (`read_records`), or `None` where its reading fails.
#[cfg(all(unix, feature = "parallel"))]
fn read_part(
    file: &File,
    bytes: Range<u64>,
    config: ReaderConfig,
    failed: &std::sync::atomic::AtomicBool,
) -> Option<(Vec<usize>, u64)> {
    let mut lengths = Vec::new();
    let lines = read_records(file, bytes, config, (0, 1), false, failed, |record, _| {
        lengths.push(record.sequence.len());
        After::Next
    })
    .ok()?;
    Some((lengths, lines))
}

/// LOSAT_X_QSTREAM's second reading of part `span` of the first reading of `file`: its records read
/// again by a new reader that starts at the part's line and record (`read_records`), each with the
/// length the first reading found (`lengths`, every input record's). The last part (`last`) is read
/// up to the file's end now and only up to its records: like the one reader of
/// `x_query_stream::State::read_one`, it does not read a record added after them, and reads letters
/// added to the last record. Returns the records read and, where the file changed since the first
/// reading, the first input record it cannot give with what changed: a record of another length, a
/// reader error, a part that holds fewer records or more (the bytes after its records hold
/// another); for the part as a whole (the input record is its first), a part that no longer starts
/// at a record start (a `>` after a `\n`), whose bytes hold a `\r` or a `>?` or cannot all be read,
/// or whose lines are not the first reading's (the parts after it start at the first reading's
/// line). No NCBI counterpart: NCBI reads the query file once, with one reader; the part's reader
/// starts where that reader is (scheduling only).
#[cfg(all(unix, feature = "parallel"))]
pub(super) fn read_again(
    file: &File,
    span: &Span,
    last: bool,
    config: ReaderConfig,
    lengths: &[usize],
) -> (Vec<FastaRecord>, Option<(usize, String)>) {
    use std::os::unix::fs::FileExt;

    let whole = |what: &str| (Vec::new(), Some((span.first, what.to_string())));
    // The first part starts where the reader starts; every other at a record start.
    if span.first > 0 {
        let mut mark = [0u8; 2];
        if file.read_exact_at(&mut mark, span.bytes.start - 1).is_err() || &mark != b"\n>" {
            return whole("its part no longer starts at a record start");
        }
    }
    let end = if last {
        match file.metadata() {
            Ok(metadata) => metadata.len(),
            Err(error) => return whole(&format!("{error}")),
        }
    } else {
        span.bytes.end
    };
    let mut records: Vec<FastaRecord> = Vec::with_capacity(span.records);
    let mut changed = None;
    let read = read_records(
        file,
        span.bytes.start..end,
        config,
        (span.line, span.first as i64 + 1),
        false,
        &std::sync::atomic::AtomicBool::new(false),
        |record, _| {
            let index = span.first + records.len();
            if records.len() == span.records {
                changed = Some((index, "the part before it holds another record".to_string()));
                return After::Refuse;
            }
            if record.sequence.len() != lengths[index] {
                changed = Some((
                    index,
                    format!(
                        "{} letters where the first reading read {}",
                        record.sequence.len(),
                        lengths[index]
                    ),
                ));
                return After::Refuse;
            }
            records.push(record);
            if last && records.len() == span.records {
                After::Stop
            } else {
                After::Next
            }
        },
    );
    let index = span.first + records.len();
    match read {
        Err(PartError::Bytes) => whole(
            "its part holds a carriage return or a gap line, or ends before the first reading's end",
        ),
        // Only a record that differs gives the part up here (no other part sets `failed`).
        Err(PartError::Refused) => (
            records,
            changed.or_else(|| Some((index, "its part was given up".to_string()))),
        ),
        Err(PartError::Reader(error)) => (records, Some((index, format!("{error}")))),
        Ok(_) if records.len() < span.records => {
            let what = if last {
                "the file ends before it"
            } else {
                "its part ends before it"
            };
            (records, Some((index, what.to_string())))
        }
        Ok(lines) if !last && lines != span.lines => whole(&format!(
            "its part holds {lines} lines where the first reading read {}",
            span.lines
        )),
        Ok(_) => (records, None),
    }
}

/// Bytes `at..end` of a file, read at their offsets (the file's position stays where it is), that
/// notes (`refused`) a `\r` or a `>?` in what it returns, and an I/O error or an end of the file
/// before `end` (which the reader would take for the end of its input).
#[cfg(all(unix, feature = "parallel"))]
struct PartReader<'a> {
    file: &'a File,
    at: u64,
    end: u64,
    /// The last byte returned (a part starts at the start of a line).
    last: u8,
    refused: &'a mut bool,
}

#[cfg(all(unix, feature = "parallel"))]
impl std::io::Read for PartReader<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        use std::os::unix::fs::FileExt;

        let left = usize::try_from(self.end - self.at).unwrap_or(usize::MAX);
        let want = buf.len().min(left);
        if want == 0 {
            return Ok(0);
        }
        let n = match self.file.read_at(&mut buf[..want], self.at) {
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
        self.at += n as u64;
        Ok(n)
    }
}

#[cfg(all(test, unix, feature = "parallel"))]
#[path = "x_query_scan_par_tests.rs"]
mod tests;
