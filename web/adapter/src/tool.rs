//! The runs of the command line's tools besides the searches (`losat_web2_tool`, docs/web/abi_v2.md
//! §4, §5; line NP WW3a, 2026-11-22; line RG-F3p7; line NP WW4b1, 2026-11-27; line NP WW4c,
//! 2026-11-28): `windowmasker`, `blast_formatter`, `makeblastdb`, `dustmasker`, `segmasker` and
//! `convert2blastmask`, with the command line's words, in the run's working directory (`run::enter_run_directory`), whose files they read by name as the command line reads
//! them in its directory. Standard input is the null device (empty), standard output goes to
//! stream 11 and standard error to stream 3, then the exit status goes to stream 2
//! (`{"status":N}`).
//!
//! `windowmasker`, `blast_formatter` and the masking tools write their files (`-out FILE`, `-logfile FILE`) to stream
//! 12, each as it is written, so that no output is held whole (D-WW3-6: human counts are 724 MB, a
//! genome masked as FASTA 3.1 GB): each write at its offset, as the system writes a file that a run
//! opens twice (`-out x -logfile x`, D-WW3-9). `makeblastdb` writes, reads back, renames and
//! removes the files of its database in the run's working directory through the file system, as on
//! the command line (D-WW4b-1): the host gives a tool run a working directory that it can write
//! (abi_v2.md §3), and hands the files that the run left there to its caller when the run ends.

use std::io::{self, Write};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Arc;

use losat::web_adapter::{
    blast_formatter_read_files, convert2blastmask_read_files, dustmasker_read_files,
    segmasker_read_files, windowmasker_read_files, FileReceiver, FileStream, RunFiles, PROGRAMS,
};

use crate::emit::{
    emit, StreamWriter, CHUNK, FRAME_BEGIN, FRAME_BYTES, FRAME_END, FRAME_TRUNCATE,
    STREAM_DIAGNOSTICS, STREAM_PRINTED, STREAM_RESPONSE, STREAM_TOOL_FILES,
};

/// The tools that a tool run runs: `windowmasker`, `blast_formatter`, `makeblastdb`, `dustmasker`,
/// `segmasker` and `convert2blastmask`.
pub const TOOLS: &[&str] = &[
    "windowmasker",
    "blast_formatter",
    "makeblastdb",
    "dustmasker",
    "segmasker",
    "convert2blastmask",
];

/// A tool that a tool run runs (`TOOLS`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum ToolKind {
    WindowMasker,
    BlastFormatter,
    Makeblastdb,
    DustMasker,
    SegMasker,
    Convert2BlastMask,
}

/// The tool `words[0]` and its words after the name; `Err` for a name that is not a tool that
/// LOSAT Web runs: another program of the command line (`PROGRAMS`: a search, which
/// `losat_web2_run` runs, or a tool that LOSAT Web does not run yet), or no program.
fn tool<'a>(words: &'a [&'a str]) -> Result<(ToolKind, &'a [&'a str]), String> {
    let (name, arguments) = words.split_first().unwrap_or((&"", &[]));
    match *name {
        "windowmasker" => Ok((ToolKind::WindowMasker, arguments)),
        "blast_formatter" => Ok((ToolKind::BlastFormatter, arguments)),
        "makeblastdb" => Ok((ToolKind::Makeblastdb, arguments)),
        "dustmasker" => Ok((ToolKind::DustMasker, arguments)),
        "segmasker" => Ok((ToolKind::SegMasker, arguments)),
        "convert2blastmask" => Ok((ToolKind::Convert2BlastMask, arguments)),
        other if PROGRAMS.contains(&other) => {
            Err(format!("{other} is not available in the browser (ABI v2)"))
        }
        other => Err(format!("unknown program '{other}'")),
    }
}

/// Runs the tool `words[0]` with the words after it (docs/web/abi_v2.md §4 `losat_web2_tool`) and
/// emits its streams: 11, 3 and 12 as the tool writes them, then the exit status on 2. `Err` for a
/// name that is not a tool that LOSAT Web runs; nothing is emitted then.
pub fn run(words: &[&str]) -> Result<(), String> {
    let (kind, arguments) = tool(words)?;
    let mut out = StreamWriter::new(STREAM_PRINTED);
    let files = || RunFiles::receiver(Arc::new(ToolFiles::for_tool(kind, arguments)));
    let status = match kind {
        ToolKind::WindowMasker => losat::web_adapter::run_windowmasker(
            arguments,
            &files(),
            &mut io::empty(),
            &mut out,
            &mut Diagnostics,
        ),
        ToolKind::BlastFormatter => losat::web_adapter::run_blast_formatter(
            arguments,
            &files(),
            &mut io::empty(),
            &mut out,
            &mut Diagnostics,
        ),
        ToolKind::DustMasker => losat::web_adapter::run_dustmasker(
            arguments,
            &files(),
            &mut io::empty(),
            &mut out,
            &mut Diagnostics,
        ),
        ToolKind::SegMasker => losat::web_adapter::run_segmasker(
            arguments,
            &files(),
            &mut io::empty(),
            &mut out,
            &mut Diagnostics,
        ),
        ToolKind::Convert2BlastMask => losat::web_adapter::run_convert2blastmask(
            arguments,
            &files(),
            &mut io::empty(),
            &mut out,
            &mut Diagnostics,
        ),
        // Through the file system of the run's working directory (D-WW4b-1); a `-logfile` file
        // elsewhere on stream 12 (`ToolFiles::writes`).
        ToolKind::Makeblastdb => losat::web_adapter::run_makeblastdb(
            arguments,
            &files(),
            &mut io::empty(),
            &mut out,
            &mut Diagnostics,
        ),
    };
    out.finish();
    emit(
        STREAM_RESPONSE,
        format!("{{\"status\":{status}}}").as_bytes(),
    );
    Ok(())
}

/// Checks the words of a tool run as the command line reads them before its tool reads its input
/// or opens a file (`losat_web2_validate` of a tool's argv, docs/web/abi_v2.md §4; line NP WW3a,
/// D-WW3-7 and D-WW3-10; D-WW4b-3): `Err` with the command line's error (which reads the
/// registry files of the run's working directory and writes no file), or for a name that is not a
/// tool that LOSAT Web runs. Emits nothing.
pub fn validate(words: &[&str]) -> Result<(), String> {
    let (kind, arguments) = tool(words)?;
    match kind {
        ToolKind::WindowMasker => losat::web_adapter::check_windowmasker(arguments),
        ToolKind::BlastFormatter => losat::web_adapter::check_blast_formatter(arguments),
        ToolKind::Makeblastdb => losat::web_adapter::check_makeblastdb(arguments),
        ToolKind::DustMasker => losat::web_adapter::check_dustmasker(arguments),
        ToolKind::SegMasker => losat::web_adapter::check_segmasker(arguments),
        ToolKind::Convert2BlastMask => losat::web_adapter::check_convert2blastmask(arguments),
    }
}

/// The standard error of a tool run: stream 3, each write as it comes (the tool's diagnostics, a
/// message at a time), in emits of at most `CHUNK` bytes (review of WW3a, L6).
struct Diagnostics;

impl Write for Diagnostics {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        for chunk in data.chunks(CHUNK) {
            emit(STREAM_DIAGNOSTICS, chunk);
        }
        Ok(data.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// The files of a tool run, which go to stream 12 as they are written (`FileReceiver::open`),
/// numbered in the order the run opens them, carrying the run's read input names so `-logfile`
/// and `-out` naming a missing input the CLI reads back are refused.
#[derive(Default)]
struct ToolFiles {
    next: AtomicU32,
    read: Vec<Vec<u8>>,
    out_read: Vec<Vec<u8>>,
    /// Whether the run writes the files of its working directory through the file system
    /// (makeblastdb, D-WW4b-1; `FileReceiver::writes_working_directory`): only its files elsewhere
    /// (`-logfile /tmp/mk.log`) come here.
    writes: bool,
}

impl ToolFiles {
    fn for_tool(kind: ToolKind, arguments: &[&str]) -> Self {
        let (read, out_read) = match kind {
            ToolKind::WindowMasker => windowmasker_read_files(arguments),
            ToolKind::BlastFormatter => blast_formatter_read_files(arguments),
            ToolKind::DustMasker => dustmasker_read_files(arguments),
            ToolKind::SegMasker => segmasker_read_files(arguments),
            ToolKind::Convert2BlastMask => convert2blastmask_read_files(arguments),
            // makeblastdb's files of the working directory go through the file system (`run`),
            // where the command line's refusal of a log naming an input applies.
            ToolKind::Makeblastdb => (Vec::new(), Vec::new()),
        };
        Self {
            next: AtomicU32::new(0),
            read,
            out_read,
            writes: kind == ToolKind::Makeblastdb,
        }
    }
}

impl FileReceiver for ToolFiles {
    /// A file given whole (no file of windowmasker's): its frames at once.
    fn receive(&self, name: &[u8], bytes: &[u8]) -> io::Result<()> {
        self.opened(name).write_at(0, bytes)
    }

    fn open(&self, name: &[u8]) -> Option<Box<dyn FileStream>> {
        Some(Box::new(self.opened(name)))
    }

    fn read_files(&self) -> &[Vec<u8>] {
        &self.read
    }

    fn output_read_files(&self) -> &[Vec<u8>] {
        &self.out_read
    }

    fn writes_working_directory(&self) -> bool {
        self.writes
    }
}

impl ToolFiles {
    /// The file `name`, whose beginning is emitted.
    fn opened(&self, name: &[u8]) -> FileFrames {
        let number = self.next.fetch_add(1, Ordering::Relaxed);
        frame(number, FRAME_BEGIN, &[name]);
        FileFrames {
            number,
            offset: 0,
            pending: Vec::new(),
        }
    }
}

/// The largest number of a file's bytes in one frame, after its header and offset, so that a
/// frame is one `emit` call.
const FRAME_PAYLOAD: usize = CHUNK - 16;

/// Emits one frame of `STREAM_TOOL_FILES`: the file's number, the kind, then the parts of its
/// payload.
fn frame(number: u32, kind: u32, payload: &[&[u8]]) {
    let mut bytes = Vec::with_capacity(8 + payload.iter().map(|part| part.len()).sum::<usize>());
    bytes.extend_from_slice(&number.to_le_bytes());
    bytes.extend_from_slice(&kind.to_le_bytes());
    for part in payload {
        bytes.extend_from_slice(part);
    }
    emit(STREAM_TOOL_FILES, &bytes);
}

/// A file of a tool run: its writes go out in frames of bytes at their offsets, those that follow
/// each other joined up to `FRAME_PAYLOAD` bytes; a write elsewhere or a change of length sends
/// what waits first, so that the host applies the run's writes in their order. The end goes out
/// when the run's files end (the stream is dropped, at the end of the run).
struct FileFrames {
    number: u32,
    /// The offset of the first byte of `pending`.
    offset: u64,
    pending: Vec<u8>,
}

impl FileFrames {
    /// Sends the bytes that wait.
    fn send(&mut self) {
        if !self.pending.is_empty() {
            frame(
                self.number,
                FRAME_BYTES,
                &[&self.offset.to_le_bytes(), &self.pending],
            );
            self.offset += self.pending.len() as u64;
            self.pending.clear();
        }
    }
}

impl FileStream for FileFrames {
    fn write_at(&mut self, offset: u64, bytes: &[u8]) -> io::Result<()> {
        if offset != self.offset + self.pending.len() as u64 {
            self.send();
            self.offset = offset;
        }
        let mut rest = bytes;
        while !rest.is_empty() {
            let take = rest.len().min(FRAME_PAYLOAD - self.pending.len());
            self.pending.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if self.pending.len() == FRAME_PAYLOAD {
                self.send();
            }
        }
        Ok(())
    }

    fn set_len(&mut self, length: u64) -> io::Result<()> {
        self.send();
        frame(self.number, FRAME_TRUNCATE, &[&length.to_le_bytes()]);
        Ok(())
    }
}

impl Drop for FileFrames {
    fn drop(&mut self) {
        self.send();
        frame(self.number, FRAME_END, &[]);
    }
}

#[cfg(test)]
#[path = "tool_tests.rs"]
mod tests;
