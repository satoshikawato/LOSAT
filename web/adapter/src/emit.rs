//! Output streams to the host (docs/web/abi_v2.md §3, §5): `losat_host.emit(stream,
//! ptr, len)`. The engine formats on the thread that called the export (the caller
//! runs the search pool's slot zero, `LOSAT/src/utils/threading.rs`), so every call
//! reaches the host's own thread.

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
