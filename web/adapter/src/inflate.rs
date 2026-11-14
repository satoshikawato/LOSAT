//! The inflation of a compressed source file (docs/web/abi_v2.md §4, §5, §9): the application
//! gives a file named `.gz`, `.bz2` or `.zst` in pieces, and the decoder that the command line
//! reads such a file with (`PushInflater`, the push form of `Compression::decoder`) inflates it.
//! A file then gives LOSAT Web the command line's bytes, or fails with the command line's
//! sentence without its `-query: ` (swebok-improve SW8 DECISIONS 44). The name rule is the
//! engine's (`compression_of_name`), so the application holds no copy of it.

use std::path::Path;

use losat::web_adapter::{compression_of_name, Compression, PushInflater};

use crate::emit::{emit, STREAM_INFLATED};

/// The name of the inflated file when nothing is left of the name without its suffix.
pub const EMPTY_PLAIN_NAME: &str = "input";

/// The inflation of one file, from `begin` to `end`.
pub struct Inflation {
    name: String,
    compression: Compression,
    inflater: PushInflater,
}

impl Inflation {
    /// The inflation of the file named `name`, with its *inflate* response (§9), or `None` when
    /// the name is plain (the command line reads such a file as it is).
    pub fn begin(name: &str) -> Option<(Self, String)> {
        let compression = compression_of_name(Path::new(name))?;
        // The suffix is ASCII and matched byte for byte (in any case), so the cut is on a
        // character boundary.
        let plain = &name[..name.len() - compression.suffix().len()];
        let plain = if plain.is_empty() {
            EMPTY_PLAIN_NAME
        } else {
            plain
        };
        let mut response = String::from("{\"compression\":");
        crate::json::string(&mut response, compression.format());
        response.push_str(",\"plain_name\":");
        crate::json::string(&mut response, plain);
        response.push('}');
        let inflation = Self {
            name: name.to_string(),
            compression,
            inflater: PushInflater::new(compression),
        };
        Some((inflation, response))
    }

    /// Inflates the next bytes of the file and emits on stream 5 every inflated byte that the
    /// decoder has released; the error is the command line's sentence.
    pub fn chunk(&mut self, bytes: &[u8]) -> Result<(), String> {
        let result = self.inflater.push(bytes, &mut emit_inflated);
        result.map_err(|error| self.compression.not_data(&self.name, &error.to_string()))
    }

    /// The end of the file: emits the rest of the inflated bytes, or fails with the sentence of
    /// data cut short or of an empty file.
    pub fn end(self) -> Result<(), String> {
        let Self {
            name,
            compression,
            inflater,
        } = self;
        let result = inflater.finish(&mut emit_inflated);
        result.map_err(|error| compression.not_data(&name, &error.to_string()))
    }
}

fn emit_inflated(bytes: &[u8]) {
    emit(STREAM_INFLATED, bytes);
}
