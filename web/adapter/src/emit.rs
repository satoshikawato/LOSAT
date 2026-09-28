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
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.buffer.extend_from_slice(data);
        self.position.fetch_add(data.len() as u64, Ordering::SeqCst);
        if self.buffer.len() >= CHUNK {
            let full = self.buffer.len() / CHUNK * CHUNK;
            emit(self.stream, &self.buffer[..full]);
            self.buffer.drain(..full);
        }
        Ok(data.len())
    }

    /// Positions are published on every write, so a flush has nothing to add; the
    /// engine flushes at every HSP boundary when an observer is present.
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
