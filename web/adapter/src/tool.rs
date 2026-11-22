//! The runs of the command line's tools besides the searches (`losat_web2_tool`, docs/web/abi_v2.md
//! §4, §5; line NP WW3a, 2026-11-22): `windowmasker`, with the command line's words, in the run's
//! working directory (`run::enter_run_directory`), whose files it reads by name as the command line
//! reads them in its directory. Its standard input is the null device (empty), its standard output
//! goes to stream 11 and its standard error to stream 3, and the files that it writes
//! (`-out FILE`, `-logfile FILE`) to stream 12, each as it is written, so that no output is held
//! whole (D-WW3-6: human counts are 724 MB, a genome masked as FASTA 3.1 GB); then its exit status
//! goes to stream 2 (`{"status":N}`).

use std::io::{self, Write};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Arc;

use losat::web_adapter::{FileReceiver, RunFiles};

use crate::emit::{
    emit, StreamWriter, CHUNK, FRAME_BEGIN, FRAME_BYTES, FRAME_END, STREAM_DIAGNOSTICS,
    STREAM_PRINTED, STREAM_RESPONSE, STREAM_TOOL_FILES,
};

/// The tools that a tool run runs: `windowmasker`.
pub const TOOLS: &[&str] = &["windowmasker"];

/// The command line's other programs, which a tool run does not run: the searches (`losat_web2_run`)
/// and the tools that LOSAT Web does not run yet.
const OTHER_PROGRAMS: &[&str] = &[
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "psiblast",
    "rpsblast",
    "rpstblastn",
    "makeblastdb",
    "blastdbcmd",
    "dustmasker",
    "segmasker",
    "convert2blastmask",
    "makembindex",
    "rpsbproc",
];

/// The tool `words[0]` and its words after the name; `Err` for a name that is not a tool that
/// LOSAT Web runs.
fn tool<'a>(words: &'a [&'a str]) -> Result<&'a [&'a str], String> {
    let (name, arguments) = words.split_first().unwrap_or((&"", &[]));
    match *name {
        "windowmasker" => Ok(arguments),
        other if OTHER_PROGRAMS.contains(&other) => {
            Err(format!("{other} is not available in the browser (ABI v2)"))
        }
        other => Err(format!("unknown program '{other}'")),
    }
}

/// Runs the tool `words[0]` with the words after it (docs/web/abi_v2.md §4 `losat_web2_tool`) and
/// emits its streams: 11, 3 and 12 as the tool writes them, then the exit status on 2. `Err` for a
/// name that is not a tool that LOSAT Web runs; nothing is emitted then.
pub fn run(words: &[&str]) -> Result<(), String> {
    let arguments = tool(words)?;
    let mut out = StreamWriter::new(STREAM_PRINTED);
    let mut err = StreamWriter::new(STREAM_DIAGNOSTICS);
    let status = losat::web_adapter::run_windowmasker(
        arguments,
        &RunFiles::Receiver(Arc::new(ToolFiles::default())),
        &mut io::empty(),
        &mut out,
        &mut err,
    );
    out.finish();
    err.finish();
    emit(
        STREAM_RESPONSE,
        format!("{{\"status\":{status}}}").as_bytes(),
    );
    Ok(())
}

/// Checks the words of a tool run as the command line reads them before its tool opens a file
/// (`losat_web2_validate` of a tool's argv, docs/web/abi_v2.md §4; line NP WW3a, D-WW3-7,
/// 2026-11-22): `Err` with the command line's error (`web_adapter::check_windowmasker`), or for a
/// name that is not a tool that LOSAT Web runs. Emits nothing.
pub fn validate(words: &[&str]) -> Result<(), String> {
    losat::web_adapter::check_windowmasker(tool(words)?)
}

/// The files of a tool run, which go to stream 12 as they are written (`FileReceiver::open`),
/// numbered in the order the run opens them.
#[derive(Default)]
struct ToolFiles {
    next: AtomicU32,
}

impl FileReceiver for ToolFiles {
    /// A file given whole (no file of windowmasker's): its frames at once.
    fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()> {
        let mut file = self.opened(name);
        file.write_all(bytes)
    }

    fn open(&self, name: &[u8]) -> Option<Box<dyn Write + Send>> {
        Some(Box::new(self.opened(name)))
    }
}

impl ToolFiles {
    /// The file `name`, whose beginning is emitted.
    fn opened(&self, name: &[u8]) -> FileFrames {
        let number = self.next.fetch_add(1, Ordering::Relaxed);
        frame(number, FRAME_BEGIN, name);
        FileFrames {
            number,
            buffer: Vec::new(),
        }
    }
}

/// The largest payload of a frame of bytes, so that a frame is one `emit` call.
const FRAME_PAYLOAD: usize = CHUNK - 8;

/// Emits one frame of `STREAM_TOOL_FILES`: the file's number, the kind, then `payload`.
fn frame(number: u32, kind: u32, payload: &[u8]) {
    let mut bytes = Vec::with_capacity(8 + payload.len());
    bytes.extend_from_slice(&number.to_le_bytes());
    bytes.extend_from_slice(&kind.to_le_bytes());
    bytes.extend_from_slice(payload);
    emit(STREAM_TOOL_FILES, &bytes);
}

/// An open file of a tool run: its bytes go out in frames of `FRAME_PAYLOAD` bytes as they fill,
/// the rest and the end when the run closes the file (drops it).
struct FileFrames {
    number: u32,
    buffer: Vec<u8>,
}

impl Write for FileFrames {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let mut rest = data;
        while !rest.is_empty() {
            let take = rest.len().min(FRAME_PAYLOAD - self.buffer.len());
            self.buffer.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if self.buffer.len() == FRAME_PAYLOAD {
                frame(self.number, FRAME_BYTES, &self.buffer);
                self.buffer.clear();
            }
        }
        Ok(data.len())
    }

    /// The bytes stay until a frame fills or the file is closed: a flush adds nothing that the
    /// host could see before the end.
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Drop for FileFrames {
    fn drop(&mut self) {
        if !self.buffer.is_empty() {
            frame(self.number, FRAME_BYTES, &self.buffer);
        }
        frame(self.number, FRAME_END, &[]);
    }
}

#[cfg(test)]
#[path = "tool_tests.rs"]
mod tests;
