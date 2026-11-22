//! Output streams to the host (docs/web/abi_v2.md §3, §5): the `losat_host.emit` import takes a
//! stream, a pointer and a length. The engine formats on the thread that called the export (the
//! caller runs the search pool's slot zero, `LOSAT/src/utils/threading.rs`), so every call reaches
//! the host's own thread.

use std::io::{self, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

/// The largest chunk of one `emit` call.
pub const CHUNK: usize = 1 << 20;

pub const STREAM_HITS: u32 = 1;
pub const STREAM_RESPONSE: u32 = 2;
pub const STREAM_DIAGNOSTICS: u32 = 3;
/// The output of the argv's own `-outfmt`: a format besides 0, 6 and 7, or 6 or 7 with fields
/// (docs/web/abi_v2.md §5, §7).
pub const STREAM_CUSTOM: u32 = 4;
/// The inflated bytes of a compressed source file (`inflate_chunk`, `inflate_end`;
/// docs/web/abi_v2.md §4, §5).
pub const STREAM_INFLATED: u32 = 5;
/// The files of a run besides its outputs (the file of `-logfile`): a ZIP archive, emitted after
/// the run, also when it failed (docs/web/abi_v2.md §5).
pub const STREAM_RUN_FILES: u32 = 8;
/// The signal that an output stream of a run failed (it failed alone, or its search failed in a
/// run of two searches), with where its error is: the host discards what it received of it
/// (`discard`; docs/web/abi_v2.md §5).
pub const STREAM_DISCARD: u32 = 9;
/// The signal of where the messages of a run go (`messages`; docs/web/abi_v2.md §5): the warnings
/// that the command line writes to its standard error without `-logfile`. Line WQA's sweep, D7,
/// 2026-11-20.
pub const STREAM_MESSAGES: u32 = 10;
/// The text that the command line prints to its standard output in place of a run, exit status 0
/// (docs/web/abi_v2.md §5): its help (`-h`, `-help`, `-help-full`, `--help`), the XML of
/// `-xmlhelp`, or a version (`-version`, `-version-full`, `-version-full-xml`,
/// `-version-full-json`). Line WQA's D14, 2026-11-20.
pub const STREAM_PRINTED: u32 = 11;
/// The files that a tool run writes (`-out FILE`, `-logfile FILE`), each as it is written
/// (docs/web/abi_v2.md §5, version 13; line NP WW3a, D-WW3-6, 2026-11-22): every `emit` is one
/// frame of 8 bytes, the file's number in the run (from 0, in the order the run opens its files)
/// and the frame's kind (`FRAME_BEGIN`, `FRAME_BYTES`, `FRAME_END`), each an unsigned 32-bit
/// little-endian number, then its payload: the file's name (UTF-8) for a beginning, bytes of the
/// file for bytes (at most `CHUNK` - 8 a frame), none for an end. A file's frames come in order:
/// its beginning, its bytes, its end, which comes when the run closes the file; the frames of two
/// open files may alternate.
pub const STREAM_TOOL_FILES: u32 = 12;
/// The frame kinds of `STREAM_TOOL_FILES`.
pub const FRAME_BEGIN: u32 = 0;
pub const FRAME_BYTES: u32 = 1;
pub const FRAME_END: u32 = 2;
/// The place "nowhere" in the signals of `STREAM_DISCARD` and `STREAM_MESSAGES` (`place`): the
/// messages of TBLASTX's `-logfile ""` or `-logfile /dev/null`, which the run writes nowhere.
pub const PLACE_NOWHERE: u32 = u32::MAX;
/// The place "the error that ended the run" in the signals of `STREAM_DISCARD`: an output that a run
/// which failed does not keep, because the command line's run of the argv does not write it (`discard`;
/// docs/web/abi_v2.md §5, version 7; line WQA's sweep, D15, 2026-11-20).
pub const PLACE_RUN: u32 = u32::MAX - 1;

#[cfg(target_arch = "wasm32")]
#[link(wasm_import_module = "losat_host")]
extern "C" {
    #[link_name = "emit"]
    fn host_emit(stream: u32, ptr: *const u8, len: u32);
}

#[cfg(not(target_arch = "wasm32"))]
thread_local! {
    /// Native builds (the adapter's own tests) collect the emitted bytes here.
    pub static EMITTED: std::cell::RefCell<Vec<(u32, Vec<u8>)>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Sends `bytes` to the host in chunks of at most `CHUNK` bytes.
pub fn emit(stream: u32, bytes: &[u8]) {
    for chunk in bytes.chunks(CHUNK) {
        #[cfg(target_arch = "wasm32")]
        // SAFETY: the host copies the bytes during the call (docs/web/abi_v2.md §5).
        unsafe {
            host_emit(stream, chunk.as_ptr(), chunk.len() as u32)
        };
        #[cfg(not(target_arch = "wasm32"))]
        EMITTED.with(|emitted| emitted.borrow_mut().push((stream, chunk.to_vec())));
    }
}

/// Tells the host that `stream` failed in this run, that it discards what it received of it, and
/// where the error of the stream is (`error`, `place`): one `emit` on `STREAM_DISCARD` of 8
/// bytes, the stream's number and then the place, each an unsigned 32-bit little-endian number,
/// for each output stream that failed, also where no chunk of it was sent. The run sends nothing
/// more of that stream (docs/web/abi_v2.md §5).
pub fn discard(stream: u32, error: u32) {
    let mut signal = [0; 8];
    signal[..4].copy_from_slice(&stream.to_le_bytes());
    signal[4..].copy_from_slice(&error.to_le_bytes());
    emit(STREAM_DISCARD, &signal);
}

/// Tells the host where the messages of this run go (`place`, from LOSAT's start of the run,
/// `CallerStart::messages_place`): one `emit` on `STREAM_MESSAGES` of 4 bytes, the place as an
/// unsigned 32-bit little-endian number, once per run whose start did not end it, before its
/// outputs (docs/web/abi_v2.md §5; line WQA's sweep, D7, 2026-11-20).
pub fn messages(place: u32) {
    emit(STREAM_MESSAGES, &place.to_le_bytes());
}

/// A place of text that a run writes (`losat::web_adapter::TextPlace`), as the signals of
/// `STREAM_DISCARD` and `STREAM_MESSAGES` give it: the stream that holds it (stream 3, or stream 8
/// for the file of `-logfile`), or `PLACE_NOWHERE`.
pub fn place(place: losat::web_adapter::TextPlace) -> u32 {
    match place {
        losat::web_adapter::TextPlace::StandardError => STREAM_DIAGNOSTICS,
        losat::web_adapter::TextPlace::LogFile => STREAM_RUN_FILES,
        losat::web_adapter::TextPlace::Nowhere => PLACE_NOWHERE,
    }
}

/// A stream that the engine writes one output to. It sends whole chunks as they fill
/// and publishes how many bytes were written, so that an observer can take byte
/// positions at any moment without a flush.
pub struct StreamWriter {
    stream: u32,
    buffer: Vec<u8>,
    position: Arc<AtomicU64>,
}

impl StreamWriter {
    pub fn new(stream: u32) -> Self {
        Self {
            stream,
            buffer: Vec::new(),
            position: Arc::new(AtomicU64::new(0)),
        }
    }

    /// The number of bytes written so far, readable while the engine writes.
    pub fn position(&self) -> Arc<AtomicU64> {
        self.position.clone()
    }

    /// Sends what is left.
    pub fn finish(&mut self) {
        emit(self.stream, &self.buffer);
        self.buffer.clear();
    }
}

impl Write for StreamWriter {
    /// Sends whole chunks as they fill. A write of a chunk or more is not copied past the chunk
    /// that the buffer completes: its whole chunks are sent from `data` (a side file of outfmt 13
    /// or 14 comes in one write). The chunks are the same as with every byte buffered.
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.position.fetch_add(data.len() as u64, Ordering::SeqCst);
        let mut rest = data;
        if !self.buffer.is_empty() {
            let take = rest.len().min(CHUNK - self.buffer.len());
            self.buffer.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if self.buffer.len() < CHUNK {
                return Ok(data.len());
            }
            emit(self.stream, &self.buffer);
            self.buffer.clear();
        }
        let full = rest.len() / CHUNK * CHUNK;
        emit(self.stream, &rest[..full]);
        self.buffer.extend_from_slice(&rest[full..]);
        Ok(data.len())
    }

    /// Positions are published on every write, so a flush has nothing to add; the
    /// engine flushes at every HSP boundary when an observer is present.
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
