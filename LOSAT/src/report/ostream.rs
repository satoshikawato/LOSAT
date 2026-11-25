//! The report's output stream as NCBI's writes it: the command line's report stream
//! (`ReportStream`), behind NCBI's buffer of the stream (`cli::NcbiBuffer`, the one model of it,
//! `cli/stdio_buffer.rs`), the stream that NCBI's formatter writes to (`NcbiOstream`) and the tie of
//! standard error to standard output (`tie_standard_output`).

use crate::cli::{NcbiBuffer, NcbiStdio, NCBI_BUFSIZ};

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:116-119
/// NCBI's region of the report where an insertion is made: what NCBI's code does where the stream's
/// exception (the `badbit` of its `exceptions` mask, which `CBlastFormat` sets) is thrown. The
/// formatter declares it (`ReportWrite::site`) before it writes there; a report starts in `Report`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Site {
    /// The exception reaches the application's `CATCH_ALL`: "BLAST failed to write output", exit
    /// 6 (`blast_app_util.hpp:252-255`).
    Report,
    /// The alignments of outfmt 0-4 (`CDisplaySeqalign::DisplaySeqalign`), from the flush of
    /// the `Alignments` heading's `endl` to the flush at their end: the vestigial `out2`, a
    /// `CObjectOStream` on the report's stream (`showalign.cpp:1931`), flushes the stream again
    /// in its destructor and posts 802.4 (`objostr.cpp:317-324`), then as `Report`. On standard
    /// output the post's `cerr` flushes the failed `cout` first (the tie, `ncbidiag.cpp:6476-6493`)
    /// and that throws inside `~CNcbiDiag`: `std::terminate`.
    Alignments,
    /// A `CBlastTabularInfo` lives (outfmt 6, 7, 10; `blast_format.cpp:773-782`, `:2234-2235`):
    /// its destructor flushes the failed stream again (`tabular.cpp:160-163`), the exception
    /// leaves an implicitly `noexcept` destructor: `std::terminate`.
    Tabular,
    /// SAM (outfmt 17): `CSAM_Formatter` holds the report and writes it in its destructor
    /// (`sam_formatter.cpp:313-316`), a member of `CBlastFormat` destroyed after the body of
    /// `~CBlastFormat` resets the stream's exceptions (`blast_format.hpp:415`): a failed write is
    /// not reported.
    Silent,
    /// One query's `Iteration` of BLAST XML (outfmt 5), from its first insertion to its last
    /// (`s_BlastXMLAddIteration`, inside the `try` of `BlastXML_FormatReport`,
    /// blastxml_format.cpp:745-764): the exception is caught there, `Failed
    /// s_BlastXMLAddIteration Q(0/1` (sic) is posted, and the run goes on with the next query,
    /// whose first insertion into the bad stream throws again. The error of the failed insertion
    /// is `ReportWriteFailure::caught` (`is_caught_report_write_failure`): the query's report stops
    /// there (no later `Hit` is made). On standard output the post's tie throws
    /// (`std::terminate`).
    XmlIteration,
    /// One object that an object stream of NCBI's writes on the report's stream (`m_Outfile` with
    /// `MSerial_AsnText` and the object, `WriteObject`, serialobject.cpp:815-827: the Seq-annot of
    /// a query in outfmt 8, 9 and 12): its buffer (`COStreamBuffer`) inserts the object in pieces
    /// (`FlushBuffer`, at the end `EndOfWrite`), and `~CObjectOStream`'s `Close` flushes the stream
    /// (`ReportWrite::ncbi_flush`). A failed insertion throws `ios_base::failure` through the
    /// object's frames (they catch `CException`s only, objstack.hpp:201-214); `Close`, which clears
    /// the stream's state for its flush (`STemporarilyClearStreamState`), fails again on a file
    /// (whose buffer keeps the bytes it could not write, as NCBI BLAST+ 2.17.0's `-out` runs show)
    /// and posts 802.4, and the exception ends the run at `CATCH_ALL`. Where only `Close`'s flush
    /// fails, its destructor's catch posts 802.4 and the run goes on with a bad stream, whose next
    /// insertion throws. On standard output each post's tie throws (`std::terminate`). `whole`:
    /// LOSAT inserts the object in one insertion, and a failure of its device write is decided by
    /// NCBI's pieces (`ReportWrite::resolve_object`, which the writer calls after an
    /// `is_undecided_report_write_failure`); else it inserts NCBI's pieces.
    Serial { whole: bool },
    /// One query's report of outfmt 16 or 15 (`json`), written by an object stream of NCBI's on
    /// the report's stream inside the `try` of `BlastXML2_FormatReport` / `BlastJSON_FormatReport`
    /// (blastxml2_format.cpp:750-763, :877-890): as `Serial` (LOSAT inserts the report whole), but
    /// a failed insertion is caught there after `Close`'s 802.4: `XML format failed` (`JSON format
    /// failed`) is posted and the run goes on (`ReportWriteFailure::caught`); where only `Close`'s
    /// flush fails, the 802.4 alone, and the report ends well.
    Xml2Report { json: bool },
}

impl Site {
    /// The region of a run's report in format `number` (`-outfmt`'s): the tabular reports'
    /// insertions are all made while a `CBlastTabularInfo` lives (outfmt 6, 7, 10 and 20; a
    /// tabular report has no prolog), SAM's after `~CBlastFormat`; every other format starts in
    /// `Report` and declares its other regions (`ReportWrite::site`).
    pub fn of_format(number: i32) -> Site {
        use crate::blastinput::app::ReportFormat;
        match ReportFormat::from_number(number) {
            Some(
                ReportFormat::Tabular
                | ReportFormat::TabularWithComments
                | ReportFormat::CommaSeparatedValues
                | ReportFormat::CommaSeparatedValuesWithHeader,
            ) => Site::Tabular,
            Some(ReportFormat::Sam) => Site::Silent,
            _ => Site::Report,
        }
    }
}

/// What NCBI does when the report's stream failed (`ReportStream::finish`).
#[derive(Clone, Debug, PartialEq, Eq)]
enum Failure {
    /// The exception reached `CATCH_ALL` after the texts `posts` (each ending in a line feed).
    Thrown { posts: String },
    /// `std::terminate`: the C++ runtime's text, then an abort.
    Terminate,
}

/// The text that NCBI BLAST+ 2.17.0 (its Linux build) writes to standard error when its report's
/// stream throws where nothing catches it (a tabular report's destructor flush), before it
/// aborts (exit status 134). Observed, not read from a source: the runs of the device rows of
/// `tests/db_search_manifest/*.device.tsv` (NCBI's standard error frozen as `<row_id>.err`).
pub const TERMINATE_TEXT: &str = "terminate called after throwing an instance of \
'std::__ios_failure'\n  what():  basic_ios::clear: iostream error\n";

// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:317-324
/// The post of `~CObjectOStream` whose `Close` flushed the failed report stream of `program`
/// (`[program]` is the application's name in NCBI's posts).
pub fn cannot_close_post(program: &str) -> String {
    format!(
        "Error: (802.4) [{program}] [Cannot close serializing output stream] Exception: \
         basic_ios::clear: iostream error\n"
    )
}

/// The error of an insertion into (or a flush of) the report's stream while NCBI's stream throws
/// (`exceptions(badbit)`) and a write of its device failed: the run stops there, as NCBI's
/// exception unwinds; `ReportStream::finish` turns it into NCBI's end (`Site`). Where NCBI's
/// formatter catches the exception (`caught`: `Site::XmlIteration`, whose post the stream has
/// made), the formatter stops that part of the report and goes on.
#[derive(Debug)]
pub struct ReportWriteFailure {
    kind: FailureKind,
}

impl ReportWriteFailure {
    /// The failure of an exception that unwinds (for the tests of the programs' error maps).
    #[cfg(test)]
    pub(crate) fn thrown() -> Self {
        ReportWriteFailure {
            kind: FailureKind::Thrown,
        }
    }
}

/// What the formatter does with a `ReportWriteFailure`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FailureKind {
    /// NCBI's exception unwinds: the run stops.
    Thrown,
    /// NCBI's formatter catches it (`Site::XmlIteration`): that part of the report stops.
    Caught,
    /// The device write of a whole object (`Site::Serial { whole: true }`) failed: NCBI's end
    /// depends on the pieces of the object (`ReportWrite::resolve_object`).
    Undecided,
}

impl std::fmt::Display for ReportWriteFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("BLAST failed to write output")
    }
}

impl std::error::Error for ReportWriteFailure {}

/// Whether `error` is the report stream's failure (`ReportWriteFailure`).
pub fn is_report_write_failure(error: &std::io::Error) -> bool {
    error
        .get_ref()
        .is_some_and(|inner| inner.is::<ReportWriteFailure>())
}

/// Whether `error` is the report stream's failure in a region whose exception NCBI's formatter
/// catches (`ReportWriteFailure::caught`, `Site::XmlIteration`): the formatter stops that part
/// of the report (its post is made) and goes on.
pub fn is_caught_report_write_failure(error: &std::io::Error) -> bool {
    failure_kind(error) == Some(FailureKind::Caught)
}

/// Whether `error` is the failed device write of an object that LOSAT inserted whole
/// (`Site::Serial { whole: true }`): the writer hands NCBI's pieces of the object to
/// `ReportWrite::resolve_object`, which ends it as NCBI's does.
pub fn is_undecided_report_write_failure(error: &std::io::Error) -> bool {
    failure_kind(error) == Some(FailureKind::Undecided)
}

fn failure_kind(error: &std::io::Error) -> Option<FailureKind> {
    error
        .get_ref()
        .and_then(|inner| inner.downcast_ref::<ReportWriteFailure>())
        .map(|failure| failure.kind)
}

fn write_failure() -> std::io::Error {
    std::io::Error::other(ReportWriteFailure {
        kind: FailureKind::Thrown,
    })
}

fn caught_failure() -> std::io::Error {
    std::io::Error::other(ReportWriteFailure {
        kind: FailureKind::Caught,
    })
}

fn undecided_failure() -> std::io::Error {
    std::io::Error::other(ReportWriteFailure {
        kind: FailureKind::Undecided,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:745-764
/// The post of `BlastXML_FormatReport`'s catch of a failed insertion of an `Iteration` of
/// `program` (`std::ios_base::failure` is no `CException`: no `what()` follows the label, which
/// NCBI writes without its `)`; each call formats one query, `Q(0/1`).
fn iteration_post(program: &str) -> String {
    format!("Error: [{program}] Failed s_BlastXMLAddIteration Q(0/1\n")
}

/// The command line's report stream (`-out` or standard output) as NCBI's writes it: NCBI's
/// buffer of the stream (`NcbiBuffer`: the buffers of the `-out` file and of standard output that
/// NCBI BLAST+ 2.17.0's Linux build is observed to keep, line-buffered on a terminal) before the
/// file or device, so that a write reaches the
/// device where NCBI's buffer writes (each `write` is one insertion of NCBI's, `insert`; a flush
/// of NCBI's is `ReportWrite::ncbi_flush`, LOSAT's flushes hold), and NCBI's stream state:
/// - a failed write of the device sets `badbit`; while the formatter lives (`CBlastFormat` sets
///   `exceptions(badbit)` and `~CBlastFormat` resets it, `finish`), that insertion and every later
///   one fail (`ReportWriteFailure`: an insertion into a bad stream throws, as NCBI's runs show),
///   and what NCBI does then is the region's (`Site`, `finish`);
/// - after `~CBlastFormat` a bad stream drops the insertions, and a failed write is not reported
///   (the exit flushes C's `stdout` without a check).
///
/// `cerr` is tied to `cout`: a write to standard error writes what standard output holds first
/// (`tie_standard_output`, `CliStderr`); where that fails while the formatter lives, the post's
/// destructor throws (`std::terminate`). The `-out` file is not tied: a text on standard error
/// does not write it.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:317-323
pub struct ReportStream(std::sync::Arc<std::sync::Mutex<Core>>);

/// What `ReportStream` holds.
struct Core {
    device: Box<dyn std::io::Write + Send>,
    /// The bytes inserted that the device has not taken: the first `out` of them NCBI's buffer has
    /// written (its writes `pending`, which LOSAT makes later, together), the others it holds.
    held: Vec<u8>,
    /// How many of the bytes `held` NCBI's buffer has written.
    out: usize,
    /// NCBI's writes of the bytes `held[..out]`, in order.
    pending: Vec<NcbiWrite>,
    /// NCBI's buffer of the stream (`None`: none, each insertion written at once). Its count of
    /// the bytes it holds is `held.len() - out` (`model`), which the insertions that it takes
    /// without a write do not keep in it.
    buffer: Option<NcbiBuffer>,
    /// An insertion of fewer bytes goes into NCBI's buffer without a write (`NcbiBuffer::room`);
    /// 0 where every insertion takes `insert_slow` (no buffer, a bad stream, a terminal).
    room: usize,
    /// `badbit`: a write of the device failed.
    bad: bool,
    /// `exceptions(badbit)`: the formatter lives.
    throws: bool,
    /// The region of the insertions (`Site`).
    site: Site,
    /// NCBI's end after the failure (`finish`).
    failure: Option<Failure>,
    /// Whether the stream is `cout` (tied to `cerr`).
    standard_output: bool,
    /// The application's name in NCBI's posts.
    program: &'static str,
    /// Whether the failure of the caught region being written was posted (`Site::XmlIteration`:
    /// one post per `Iteration`).
    region_posted: bool,
    /// Where NCBI's posts during the report go (`post`): the process's standard error
    /// (`PROCESS_CERR`), or a test's writer and its `cerr`.
    diagnostics: Option<(Box<dyn std::io::Write + Send>, Cerr)>,
    /// The bytes inserted so far (the stream's position).
    inserted: u64,
    /// The bytes that the device took.
    landed: u64,
    /// The failed device write of a whole object (`Site::Serial { whole: true }`), until
    /// `resolve_object`: NCBI's buffer before the object and the stream's position there.
    undecided: Option<(NcbiBuffer, u64)>,
    /// NCBI's buffer and the stream's position where the object being written started (the
    /// declaration of its region, `set_site`): LOSAT may insert a whole object in several writes
    /// (the archive's), all of which before the failed one reached the device.
    object_start: Option<(NcbiBuffer, u64)>,
    /// Every write of NCBI's is made where NCBI makes it (the tests' reference of `pending`).
    #[cfg(test)]
    immediate: bool,
}

/// One write of NCBI's buffer that LOSAT has not made yet (`Core::pending`): where it ends in the
/// stream, and what decides NCBI's end where it fails (`Core::failed_write`).
#[derive(Clone, Copy, Debug)]
struct NcbiWrite {
    /// The stream's position at its end.
    end: u64,
    /// The region of the insertion or flush that writes.
    site: Site,
    /// A flush's write (else an insertion's).
    flushing: bool,
    /// The formatter lived (`exceptions(badbit)`).
    throws: bool,
    /// The flush of `cerr`'s tie before a post (`Core::tie`): its failure terminates.
    tie: bool,
    /// The insertion of an object that LOSAT inserts whole (`Site::Serial { whole: true }`,
    /// `Site::Xml2Report`): NCBI's buffer and the stream's position where the object started.
    object: Option<(NcbiBuffer, u64)>,
}

/// How many bytes of NCBI's writes LOSAT keeps before it makes them (`Core::pending`): as many as
/// the buffer before the device held before NCBI's buffers were modelled.
const PENDING_BYTES: usize = 8192;

/// An insertion of which NCBI's buffer writes this many bytes or more is written at once, not
/// copied (with what is pending before it).
const DIRECT_BYTES: usize = 4096;

impl Site {
    /// Whether NCBI's writes in this region may be made later, together, up to the next point
    /// where something else is written (`Core::sync`): where its failure ends the run or is not
    /// reported, and BLAST XML's `Iteration`, whose end makes them (its failure is caught there).
    /// An object that LOSAT inserts whole (`Serial`, `Xml2Report`) is written where NCBI writes,
    /// so that its pieces decide its failure at its insertion.
    fn defers(self) -> bool {
        matches!(
            self,
            Site::Report | Site::Alignments | Site::Tabular | Site::Silent | Site::XmlIteration
        )
    }

    /// Whether a failure in this region ends the run (`Failure`, no catch): NCBI's writes stay
    /// pending from one such region to another.
    fn ends_run(self) -> bool {
        matches!(self, Site::Report | Site::Alignments | Site::Tabular)
    }
}

/// Writes `data` to `device`, counting the bytes it takes in `landed` (a write that stops at a
/// file's size limit takes part of the bytes, the next fails).
fn write_device(
    device: &mut dyn std::io::Write,
    landed: &mut u64,
    data: &[u8],
) -> std::io::Result<()> {
    write_device2(device, landed, data, &[])
}

/// Writes `first` and then `second` to `device` (as one `writev` where the device takes it),
/// counting the bytes it takes in `landed` (a write that stops at a file's size limit takes part
/// of the bytes, the next fails).
fn write_device2(
    device: &mut dyn std::io::Write,
    landed: &mut u64,
    mut first: &[u8],
    mut second: &[u8],
) -> std::io::Result<()> {
    while !first.is_empty() || !second.is_empty() {
        let written = if first.is_empty() {
            device.write(second)
        } else if second.is_empty() {
            device.write(first)
        } else {
            device.write_vectored(&[std::io::IoSlice::new(first), std::io::IoSlice::new(second)])
        };
        match written {
            Ok(0) => return Err(std::io::ErrorKind::WriteZero.into()),
            Ok(n) => {
                *landed += n as u64;
                let from_first = n.min(first.len());
                first = &first[from_first..];
                second = &second[n - from_first..];
            }
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(())
}

/// NCBI's insertions `pieces` of an object (an object stream's buffer writes) into NCBI's buffer
/// `buffer`, from the stream's position `position`, against a device that took the stream's first
/// `landed` bytes: the first piece at which the buffer writes past them (the write that fails), by
/// its index and the stream's position at its end; `None` where every write of the buffer inside
/// the object succeeds (`Core::resolve_object`, `side_file_failure`).
fn failed_piece(
    mut buffer: NcbiBuffer,
    mut position: u64,
    pieces: &[usize],
    landed: u64,
) -> Option<(usize, u64)> {
    for (k, &piece) in pieces.iter().enumerate() {
        position += piece as u64;
        let held = buffer.held();
        // The pieces' bytes are not read: a file's buffer counts them.
        buffer.sputn(piece, &[]);
        let writes = piece > 0 && buffer.held() < held + piece;
        if writes && position - buffer.held() as u64 > landed {
            return Some((k, position));
        }
    }
    None
}

impl Core {
    /// Whether NCBI's writes of the region being written may be made later (`Site::defers`).
    fn defers(&self) -> bool {
        #[cfg(test)]
        if self.immediate {
            return false;
        }
        self.site.defers()
    }

    /// One insertion (`NcbiBuffer::sputn`). One that NCBI's buffer takes without a write is held
    /// here; the others take `insert_slow`.
    #[inline]
    fn insert(&mut self, buf: &[u8]) -> std::io::Result<()> {
        let n = buf.len();
        if n < self.room {
            self.room -= n;
            self.inserted += n as u64;
            self.held.extend_from_slice(buf);
            return Ok(());
        }
        self.insert_slow(buf)
    }

    /// NCBI's buffer as the bytes it holds make it.
    fn model(&mut self) -> Option<NcbiBuffer> {
        let held = self.held.len() - self.out;
        self.buffer.as_mut().map(|buffer| {
            buffer.set_held(held);
            *buffer
        })
    }

    /// `room` after a change of NCBI's buffer or of the stream's state.
    fn update_room(&mut self) {
        self.room = match self.model() {
            Some(buffer) if !self.bad => buffer.room(),
            _ => 0,
        };
    }

    /// An insertion that NCBI's buffer may write: what it writes of the bytes held and `buf` is
    /// NCBI's write (`NcbiWrite`), made now in a region that does not defer it, where it is large,
    /// or where `PENDING_BYTES` are pending; it holds the rest.
    #[inline(never)]
    fn insert_slow(&mut self, buf: &[u8]) -> std::io::Result<()> {
        // An empty write is no insertion of NCBI's (a part of a report that NCBI does not write).
        if buf.is_empty() {
            return Ok(());
        }
        if self.bad {
            return self.insert_bad();
        }
        let start = self.inserted;
        self.inserted += buf.len() as u64;
        let Some(before) = self.model() else {
            // No buffer of NCBI's: the insertion is written at once.
            self.pending.push(self.write_here(false, None));
            self.out = self.held.len();
            return self.emit(buf);
        };
        let mut buffer = before;
        buffer.sputn(buf.len(), buf);
        let kept = buffer.held();
        self.buffer = Some(buffer);
        if kept == before.held() + buf.len() {
            // NCBI's buffer takes it (a terminal's insertion without a newline).
            self.held.extend_from_slice(buf);
            self.update_room();
            return Ok(());
        }
        // NCBI's buffer writes all that it held, then the insertion's first bytes.
        let object = matches!(
            self.site,
            Site::Serial { whole: true } | Site::Xml2Report { .. }
        )
        .then(|| self.object_start.unwrap_or((before, start)));
        let mut write = self.write_here(false, object);
        write.end = self.inserted - kept as u64;
        self.pending.push(write);
        let written = buf.len() - kept;
        if self.defers() && written < DIRECT_BYTES {
            self.held.extend_from_slice(buf);
            self.out = self.held.len() - kept;
            if self.out >= PENDING_BYTES {
                self.emit(&[])?;
            }
        } else {
            self.out = self.held.len();
            self.emit(&buf[..written])?;
            self.held.extend_from_slice(&buf[written..]);
        }
        self.update_room();
        Ok(())
    }

    /// An insertion into the bad stream.
    #[cold]
    #[inline(never)]
    fn insert_bad(&mut self) -> std::io::Result<()> {
        // An object written to the bad stream after `~CBlastFormat` (the error archive of
        // `PrintErrorArchive`): its buffer's first write fails without an exception mask, a
        // write fault at the buffer's position then (`write_fault`).
        if !self.throws && self.site == (Site::Serial { whole: true }) && self.undecided.is_none() {
            if let Some(buffer) = self.buffer {
                self.undecided = Some((buffer, self.inserted));
                return Err(undecided_failure());
            }
        }
        self.after_failure()
    }

    /// The insertions `report[from..at]` for each cut `(at, _)` of `cuts` (all of them
    /// `Cut::Insertion`) and `report[last..to]` (`ReportWrite::insert_cut_run`): NCBI's writes
    /// among them are made together, at the end, where the region defers them.
    fn insert_cut_run(
        &mut self,
        report: &[u8],
        from: usize,
        cuts: &[(usize, Cut)],
        to: usize,
    ) -> std::io::Result<()> {
        let defers = self.defers();
        let Some(mut buffer) = self.model().filter(|_| !self.bad && defers) else {
            let mut at = from;
            for end in cuts.iter().map(|&(end, _)| end).chain([to]) {
                let end = end.max(at);
                self.insert(&report[at..end])?;
                at = end;
            }
            return Ok(());
        };
        let start = self.inserted;
        let mut at = from;
        // The bytes of the run that NCBI's buffer has written at its last write in the run.
        let mut written = None;
        for end in cuts.iter().map(|&(end, _)| end).chain([to]) {
            let end = end.max(at);
            let n = end - at;
            if n > 0 {
                let held = buffer.held();
                // An insertion within the room is held whole (`NcbiBuffer::room`).
                buffer.sputn(n, &report[..end]);
                if buffer.held() != held + n {
                    let mut write = self.write_here(false, None);
                    write.end = start + (end - from - buffer.held()) as u64;
                    self.pending.push(write);
                    written = Some(end - from - buffer.held());
                }
            }
            at = end;
        }
        let run = &report[from..to];
        self.inserted = start + run.len() as u64;
        self.buffer = Some(buffer);
        match written {
            None => self.held.extend_from_slice(run),
            Some(written) if written >= DIRECT_BYTES => {
                self.out = self.held.len();
                self.emit(&run[..written])?;
                self.held.extend_from_slice(&run[written..]);
            }
            Some(written) => {
                self.held.extend_from_slice(run);
                self.out = self.held.len() - (run.len() - written);
                if self.out >= PENDING_BYTES {
                    self.emit(&[])?;
                }
            }
        }
        self.update_room();
        Ok(())
    }

    /// The bytes that NCBI's buffer has written (those the device took and those pending).
    #[cfg(test)]
    fn ncbi_written(&self) -> u64 {
        self.landed + self.out as u64
    }

    /// NCBI's write at the stream's position now, in the region being written.
    fn write_here(&self, flushing: bool, object: Option<(NcbiBuffer, u64)>) -> NcbiWrite {
        NcbiWrite {
            end: self.inserted,
            site: self.site,
            flushing,
            throws: self.throws,
            tie: false,
            object,
        }
    }

    /// Makes NCBI's pending writes (`held[..out]`), then `extra` (the first bytes of the insertion
    /// being made, which NCBI's buffer writes with them), to the device.
    fn emit(&mut self, extra: &[u8]) -> std::io::Result<()> {
        let Core {
            device,
            landed,
            held,
            out,
            ..
        } = self;
        let result = write_device2(&mut **device, landed, &held[..*out], extra)
            .and_then(|()| device.flush());
        match result {
            Ok(()) => {
                self.held.drain(..self.out);
                self.out = 0;
                self.pending.clear();
                Ok(())
            }
            Err(error) => self.write_failed(error),
        }
    }

    /// Makes NCBI's pending writes, before anything else is written where NCBI's would come after
    /// them (a post, another region, the run's end).
    fn sync(&mut self) -> std::io::Result<()> {
        if self.out == 0 {
            return Ok(());
        }
        self.emit(&[])
    }

    /// A write of the device that failed: NCBI's write that failed is the first that ends past
    /// the bytes that the device took, and NCBI's end is that write's (`failed_ncbi_write`).
    #[cold]
    #[inline(never)]
    fn write_failed(&mut self, error: std::io::Error) -> std::io::Result<()> {
        let landed = self.landed;
        let write = self
            .pending
            .iter()
            .find(|write| write.end > landed)
            .or(self.pending.last())
            .copied();
        self.held.clear();
        self.out = 0;
        self.pending.clear();
        self.room = 0;
        match write {
            Some(write) => self.failed_ncbi_write(write, error),
            None => {
                self.bad = true;
                Ok(())
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/connect/ncbi_socket.c:866-873 (NCBI keeps the default
    // action of SIGPIPE unless SOCK starts, which only the usage report does, after the run)
    /// NCBI's write `write` failed: the stream is bad, and NCBI's end is the region's (`Site`). A
    /// pipe whose reader closed it is such a write too (`EPIPE`: Rust's runtime ignores
    /// SIGPIPE). NCBI, whose SIGPIPE keeps its default action, is ended by it (exit status 141,
    /// no text) in every format, which LOSAT cannot be without a native call: it ends as a failed
    /// write of the format ends (approved exception 5 of `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`;
    /// line RG's T6-Q1 = A, 2026-11-17).
    fn failed_ncbi_write(
        &mut self,
        write: NcbiWrite,
        _error: std::io::Error,
    ) -> std::io::Result<()> {
        self.bad = true;
        if !write.throws {
            return Ok(());
        }
        if let (Some(object), false) = (write.object, self.standard_output) {
            self.undecided = Some(object);
            return Err(undecided_failure());
        }
        if write.tie {
            self.failure = Some(Failure::Terminate);
            return Err(write_failure());
        }
        self.thrown(write.site, write.flushing)
    }

    /// A flush of NCBI's (`flush`, `endl`, `cerr`'s tie): the buffer writes what it holds.
    fn flush(&mut self) -> std::io::Result<()> {
        if self.bad {
            return self.after_failure();
        }
        self.flush_buffer(false);
        if !self.defers() || self.out >= PENDING_BYTES {
            self.sync()?;
        }
        self.update_room();
        Ok(())
    }

    /// NCBI's buffer writes what it holds (a flush; `tie`: `cerr`'s): that write is pending. In
    /// an object that LOSAT inserts in its own pieces (`Site::Serial { whole: true }`: the archive
    /// in several writes), NCBI's pieces decide where NCBI's writes failed (`resolve_object`):
    /// the flush's failure is undecided too.
    fn flush_buffer(&mut self, tie: bool) {
        if let Some(buffer) = &mut self.buffer {
            buffer.flush();
        }
        if self.held.len() > self.out {
            let object = matches!(
                self.site,
                Site::Serial { whole: true } | Site::Xml2Report { .. }
            )
            .then_some(self.object_start)
            .flatten();
            let mut write = self.write_here(true, object);
            write.tie = tie;
            self.pending.push(write);
            self.out = self.held.len();
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-383
    /// The write fault of an object written to the bad stream after the formatter
    /// (`undecided` without the exception mask): `Close` posts 802.5 at `position`, the bytes that
    /// the object stream had put into its buffer at its first write, with an empty stack; the
    /// object is not written.
    #[cold]
    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        if self.undecided.take().is_none() {
            return Ok(());
        }
        let program = self.program;
        self.post(&format!(
            "Error: (802.5) [{program}] CObjectOStream: error at byte {position}: stack is \
             empty: write fault\n"
        ))
    }

    /// An insertion that sets the stream's `badbit` itself (`ReportWrite::badbit`): the stream
    /// throws, and NCBI's end is the application's `CATCH_ALL`; the bytes held reach the device
    /// at `finish`, as NCBI's stdio does at exit.
    #[cold]
    fn badbit(&mut self) -> std::io::Result<()> {
        self.failure.get_or_insert(Failure::Thrown {
            posts: String::new(),
        });
        Err(write_failure())
    }

    /// NCBI's end of an object whose whole insertion failed at the device
    /// (`Site::Serial { whole: true }`, `undecided`), from NCBI's pieces of it (`pieces`): from
    /// NCBI's buffer before the object, a piece that the buffer writes past the bytes the device
    /// took fails inside the object (802.4 at `Close`, then `CATCH_ALL`); where every write inside
    /// the object succeeds, only `Close`'s flush fails (802.4 posted, the run goes on).
    #[cold]
    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        let Some((buffer, position)) = self.undecided.take() else {
            return Ok(());
        };
        if failed_piece(buffer, position, pieces, self.landed).is_some() {
            if let Site::Xml2Report { json } = self.site {
                return self.xml2_caught(json);
            }
            self.failure.get_or_insert(Failure::Thrown {
                posts: cannot_close_post(self.program),
            });
            return Err(write_failure());
        }
        self.post(&cannot_close_post(self.program))
    }

    /// The exception of a bad stream at an insertion (or a flush, `flushing`) while the formatter
    /// lives, in the region `site`: what NCBI does with it.
    #[cold]
    #[inline(never)]
    fn thrown(&mut self, site: Site, flushing: bool) -> std::io::Result<()> {
        let failure = match site {
            Site::Report => Failure::Thrown {
                posts: String::new(),
            },
            // On standard output a failed write leaves nothing held (NCBI BLAST+ 2.17.0 posts no
            // 802.4 there, line RG survey T6 D3), so `out2`'s flush during the unwinding of an
            // insertion's exception, which clears the stream's state first, writes nothing and
            // succeeds: no post. Its own flush at the end of the alignments fails and posts, and
            // the post's tie throws. A file keeps the bytes it could not write: the second flush
            // fails too (the 802.4 of NCBI's `-out` runs).
            Site::Alignments if self.standard_output && !flushing => Failure::Thrown {
                posts: String::new(),
            },
            Site::Alignments if self.standard_output => Failure::Terminate,
            Site::Alignments => Failure::Thrown {
                posts: cannot_close_post(self.program),
            },
            Site::Tabular => Failure::Terminate,
            Site::Silent => {
                self.throws = false;
                return Ok(());
            }
            // The catch's post, once for the `Iteration`; the run goes on.
            Site::XmlIteration => {
                if !std::mem::replace(&mut self.region_posted, true) {
                    self.post(&iteration_post(self.program))?;
                }
                return Err(caught_failure());
            }
            Site::Serial { .. } if self.standard_output => Failure::Terminate,
            // `Close`'s flush alone failed: its destructor posts and the run goes on.
            Site::Serial { .. } if flushing => {
                self.post(&cannot_close_post(self.program))?;
                return Ok(());
            }
            Site::Serial { .. } => Failure::Thrown {
                posts: cannot_close_post(self.program),
            },
            Site::Xml2Report { .. } if self.standard_output => Failure::Terminate,
            Site::Xml2Report { .. } if flushing => {
                self.post(&cannot_close_post(self.program))?;
                return Ok(());
            }
            Site::Xml2Report { json } => return self.xml2_caught(json),
        };
        // The first failure decides NCBI's end; a later one of the same exception (an insertion
        // after it, while the run unwinds) adds nothing.
        self.failure.get_or_insert(failure);
        Err(write_failure())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:884-890
    /// The catch of a report of outfmt 16 or 15 (`json`) whose insertion failed: `Close`'s 802.4
    /// while the exception unwinds, then the catch's post; the run goes on.
    #[cold]
    fn xml2_caught(&mut self, json: bool) -> std::io::Result<()> {
        let format = if json { "JSON" } else { "XML" };
        self.post(&format!(
            "{}Error: [{}] {format} format failed\n",
            cannot_close_post(self.program),
            self.program
        ))?;
        Err(caught_failure())
    }

    /// An insertion or flush of a bad stream: it throws while the formatter lives (the region's
    /// end, `thrown`: the first insertion of a query's `Iteration` after a caught failure throws
    /// again).
    fn after_failure(&mut self) -> std::io::Result<()> {
        if self.throws {
            self.thrown(self.site, false)
        } else {
            Ok(())
        }
    }

    /// A post of NCBI's on `cerr` (`ERR_POST`) while the report is written: `cerr` is tied to
    /// `cout`, so a report on standard output writes what it holds first; where that stream
    /// failed while the formatter lives, the tie throws inside the post's destructor
    /// (`std::terminate`) and the text is not written (`tie`). Before a post on a report in a
    /// file NCBI's writes before it are made (`sync`): where one failed, NCBI's run ended there,
    /// without the post.
    /// A failed `cerr` takes no text and flushes no tie (`Cerr`).
    fn post(&mut self, text: &str) -> std::io::Result<()> {
        let cerr = match &self.diagnostics {
            Some((_, cerr)) => cerr,
            None => &PROCESS_CERR,
        };
        if text.is_empty() || cerr.failed() {
            return Ok(());
        }
        if self.standard_output {
            self.tie()?;
        } else {
            self.sync()?;
        }
        match &mut self.diagnostics {
            Some((diagnostics, cerr)) => cerr.write(&mut **diagnostics, text.as_bytes()),
            None => PROCESS_CERR.write_process(text.as_bytes()),
        }
        Ok(())
    }

    /// The region of the insertions from here on (`ReportWrite::site`): each declaration starts a
    /// region (the next query's `Iteration` posts its failure again). NCBI's writes of a region
    /// whose failure does not end the run, or before one, are made at its ends.
    /// Where those writes fail in a region whose failure NCBI's formatter catches, that region
    /// ended there (its post made) and the new one starts all the same: the caught failure is
    /// returned after the change (`is_caught_report_write_failure`).
    fn set_site(&mut self, site: Site) -> std::io::Result<()> {
        let synced = if self.site.ends_run() && site.ends_run() {
            Ok(())
        } else {
            self.sync()
        };
        if let Err(error) = synced {
            if !is_caught_report_write_failure(&error) {
                return Err(error);
            }
            self.change_site(site);
            return Err(error);
        }
        self.change_site(site);
        Ok(())
    }

    /// `set_site`'s change of the region.
    fn change_site(&mut self, site: Site) {
        self.region_posted = false;
        self.site = site;
        let inserted = self.inserted;
        self.object_start = match site {
            Site::Serial { whole: true } | Site::Xml2Report { .. } => {
                self.model().map(|buffer| (buffer, inserted))
            }
            _ => None,
        };
    }

    /// `cerr`'s tie: standard output writes what it holds before a text on standard error. A
    /// failed `cout` while the formatter lives throws there, inside the post's destructor
    /// (`std::terminate`); the text is not written.
    fn tie(&mut self) -> std::io::Result<()> {
        if self.bad {
            // A flush of a bad `cout` throws while the formatter lives (NCBI's runs terminate).
            if self.throws {
                self.failure = Some(Failure::Terminate);
                return Err(write_failure());
            }
            return Ok(());
        }
        self.flush_buffer(true);
        let tied = self.sync();
        self.update_room();
        tied
    }

    /// `~CBlastFormat` and the application's `CATCH_ALL`: the stream no longer throws and
    /// writes what it holds (a failure unreported); the run's end is NCBI's for the failure, or
    /// `result`.
    fn finish(&mut self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        let _ = self.sync();
        // A tabular report's bytes reach the device while its `CBlastTabularInfo` lives (the
        // destructor of each query's, and of the epilog's, flushes): what is held fails there.
        if self.site == Site::Tabular {
            let _ = self.flush().and_then(|()| self.sync());
        }
        self.throws = false;
        let _ = self.flush().and_then(|()| self.sync());
        match self.failure.take() {
            None => result,
            Some(Failure::Thrown { posts }) => Err(crate::cli::NativeError {
                exit: 6,
                message: [posts.as_bytes(), b"BLAST failed to write output\n"].concat(),
            }
            .into()),
            Some(Failure::Terminate) => Err(crate::cli::NativeError {
                exit: -6,
                message: TERMINATE_TEXT.as_bytes().to_vec(),
            }
            .into()),
        }
    }
}

impl Drop for Core {
    fn drop(&mut self) {
        if !self.bad {
            let _ = self.device.write_all(&self.held);
            let _ = self.device.flush();
        }
    }
}

/// The process's standard output as the command line's report device: the command line's one
/// unbuffered device of it (`cli::standard_output`), so that each write is one write of the device
/// and the bytes that a failed write left are known (`Core::write_failed`; `std::io::Stdout` keeps
/// a line's tail until its flush and counts it as taken).
fn standard_output_device() -> Box<dyn std::io::Write + Send> {
    Box::new(crate::cli::standard_output())
}

/// The run's report stream, which a write to standard error makes NCBI's writes of first
/// (`tie_standard_output`): standard output (the process's, tied to standard error) or the
/// command line's `-out` file.
static TIED_REPORT: std::sync::Mutex<Option<std::sync::Weak<std::sync::Mutex<Core>>>> =
    std::sync::Mutex::new(None);

/// The application's name in NCBI's posts of the command line's run (`ReportStream::of_cli` in
/// the command line's process, `cli::is_command_line_run`), for the posts of the streams that NCBI
/// writes apart from the report (`write_side_file`).
static RUN_PROGRAM: std::sync::Mutex<Option<&'static str>> = std::sync::Mutex::new(None);

#[inline]
fn lock<T>(mutex: &std::sync::Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

impl ReportStream {
    fn new(
        device: Box<dyn std::io::Write + Send>,
        buffer: Option<NcbiBuffer>,
        standard_output: bool,
        program: &'static str,
        site: Site,
    ) -> Self {
        // The bytes held: NCBI's buffer's and the writes pending before it (`PENDING_BYTES`).
        let held = match buffer {
            Some(_) => Vec::with_capacity(2 * NCBI_BUFSIZ),
            None => Vec::new(),
        };
        let room = buffer.map_or(0, |buffer| buffer.room());
        ReportStream(std::sync::Arc::new(std::sync::Mutex::new(Core {
            device,
            held,
            out: 0,
            pending: Vec::new(),
            buffer,
            room,
            bad: false,
            throws: true,
            site,
            failure: None,
            standard_output,
            program,
            region_posted: false,
            diagnostics: None,
            inserted: 0,
            landed: 0,
            undecided: None,
            object_start: None,
            #[cfg(test)]
            immediate: false,
        })))
    }

    /// The command line's report stream of `program` (the application's name in NCBI's posts)
    /// for the report of `-outfmt` `outfmt` (`Site::of_format`): the `-out` file `out_file`
    /// (NCBI's file buffer, `NcbiStdio::File`: 8,190 bytes held, the 8,191st written with them)
    /// or, with `None`, standard output (`NcbiStdio::cout`, tied to standard error).
    ///
    /// A run whose strategy file `strategy` LOSAT refuses gives it the stream's request point
    /// (`RequestPoint`), where NCBI's writes so far are made before the refusal.
    pub fn of_cli(
        out_file: Option<std::fs::File>,
        program: &'static str,
        outfmt: &str,
        strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    ) -> Self {
        let number = crate::blastinput::app::parse_formatting_string(outfmt)
            .map_or(0, |choice| choice.number);
        let stream = match out_file {
            Some(file) => Self::file(Box::new(file), program, Site::of_format(number)),
            None => {
                Self::standard_output(standard_output_device(), program, Site::of_format(number))
            }
        };
        if let Some(strategy) = strategy {
            strategy.attach_report(RequestPoint(std::sync::Arc::downgrade(&stream.0)));
        }
        if crate::cli::is_command_line_run() {
            *lock(&RUN_PROGRAM) = Some(program);
            *lock(&TIED_REPORT) = Some(std::sync::Arc::downgrade(&stream.0));
        }
        stream
    }

    /// The `-out` file `file` of the library's `run_to_file` (`api::local_blast::OutFile`): the
    /// command line's `-out` stream of `program` for `-outfmt` `outfmt` (`of_cli`), whose posts
    /// go to `posts` (the call's diagnostics, `api::local_blast::Posts`) where the command line
    /// writes them to its standard error.
    pub(crate) fn of_library_file(
        file: std::fs::File,
        program: &'static str,
        outfmt: &str,
        posts: Box<dyn std::io::Write + Send>,
    ) -> Self {
        let stream = Self::of_cli(Some(file), program, outfmt, None);
        lock(&stream.0).diagnostics = Some((posts, Cerr::new()));
        stream
    }

    /// The `-out` file `inner`, its report starting in `site`.
    pub fn file(inner: Box<dyn std::io::Write + Send>, program: &'static str, site: Site) -> Self {
        Self::new(
            inner,
            Some(NcbiBuffer::new(NcbiStdio::File)),
            false,
            program,
            site,
        )
    }

    /// A report kept in memory (`inner`, TBLASTN's Web ABI v1): each insertion is written at once
    /// (no stream of NCBI's to model).
    pub fn memory(inner: Box<dyn std::io::Write + Send>, program: &'static str) -> Self {
        Self::new(inner, None, false, program, Site::Report)
    }

    /// Standard output, writing to `inner` (the process's standard output), tied to standard
    /// error (`tie_standard_output`).
    ///
    /// NCBI writes the report to `cout`; when the process starts with its standard output
    /// closed (`>&-`), its first write fails (outfmt 0: "BLAST failed to write output", exit
    /// 6; outfmt 6 and 7 abort). Rust's runtime opens `/dev/null` on a closed standard
    /// descriptor before `main`, which LOSAT cannot tell from a `/dev/null` that the caller
    /// opened (read and write, as Python's `subprocess.DEVNULL`), so the report is written
    /// there and the run succeeds (approved exception 6 of
    /// `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, accepted by the maintainer in session S08b).
    pub fn standard_output(
        inner: Box<dyn std::io::Write + Send>,
        program: &'static str,
        site: Site,
    ) -> Self {
        // NCBI's Windows (MSVC) and macOS (libc++, BSD stdio) builds buffer otherwise: there
        // the exit status of a failing write is unverified (this model is their Linux builds').
        let report = Self::new(
            inner,
            Some(NcbiBuffer::new(NcbiStdio::cout())),
            true,
            program,
            site,
        );
        *lock(&TIED_REPORT) = Some(std::sync::Arc::downgrade(&report.0));
        report
    }

    /// The handle by which a formatter that writes into the stream declares NCBI's regions and
    /// flushes (`NcbiOstream::report`, `ReportWrite`).
    pub fn handle(&self) -> ReportHandle {
        ReportHandle(self.0.clone())
    }

    /// `~CBlastFormat` and the application's `CATCH_ALL`, after the run's `result`: the stream
    /// stops throwing and writes what it holds, unchecked; a failure of the stream while the
    /// formatter lived ends the run as NCBI's does (`Site`: the region's posts and "BLAST failed
    /// to write output", exit 6, or `std::terminate`'s text and an abort, `NativeError` of exit
    /// -6), whatever `result` is (the failure's error, `ReportWriteFailure`, or what the run made
    /// of it); otherwise `result`.
    pub fn finish(&mut self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        lock(&self.0).finish(result)
    }

    /// Whether a write of the stream's device failed (`badbit`), once NCBI's writes so far are
    /// made.
    pub fn failed(&self) -> bool {
        let mut core = lock(&self.0);
        let _ = core.sync();
        core.bad
    }
}

impl std::io::Write for ReportStream {
    /// One insertion of NCBI's.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        lock(&self.0).insert(buf)?;
        Ok(buf.len())
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        lock(&self.0).insert(buf)
    }

    /// LOSAT's flushes of the report hold (NCBI's are `ReportWrite::ncbi_flush`; NCBI's writes
    /// are made before a text on standard error, `tie_standard_output`).
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The command line's report stream where NCBI writes its first search strategy request
/// (`archive::strategy::check_report_stream`), for a strategy file that LOSAT refuses.
pub struct RequestPoint(std::sync::Weak<std::sync::Mutex<Core>>);

impl RequestPoint {
    /// Makes NCBI's writes so far, while the formatter lives: a failed write ends the run as
    /// NCBI's run ends there, before its request. A stream already gone has nothing to make.
    pub fn make_writes(&self) -> std::io::Result<()> {
        match self.0.upgrade() {
            Some(core) => lock(&core).sync(),
            None => Ok(()),
        }
    }
}

/// A formatter's hold on the command line's report stream (`ReportStream::handle`), by which it
/// declares NCBI's regions (`Site`) and flushes; the stream of a library caller or LOSAT Web has
/// none (`NcbiOstream::report`).
#[derive(Clone)]
pub struct ReportHandle(std::sync::Arc<std::sync::Mutex<Core>>);

impl ReportHandle {
    fn site(&self, site: Site) -> std::io::Result<()> {
        lock(&self.0).set_site(site)
    }

    #[inline]
    fn insert(&self, buf: &[u8]) -> std::io::Result<()> {
        lock(&self.0).insert(buf)
    }

    fn flush(&self) -> std::io::Result<()> {
        lock(&self.0).flush()
    }

    fn resolve_object(&self, pieces: &[usize]) -> std::io::Result<()> {
        lock(&self.0).resolve_object(pieces)
    }

    fn write_fault(&self, position: usize) -> std::io::Result<()> {
        lock(&self.0).write_fault(position)
    }

    fn badbit(&self) -> std::io::Result<()> {
        lock(&self.0).badbit()
    }

    /// The insertions from here on are made in NCBI's region `site` (`ReportWrite::site`, for a
    /// hold that keeps the regions in order with its bytes).
    pub(crate) fn declare(&self, site: Site) -> std::io::Result<()> {
        self.site(site)
    }

    /// A flush of NCBI's (`ReportWrite::ncbi_flush`, for a hold).
    pub(crate) fn ncbi_flush(&self) -> std::io::Result<()> {
        self.flush()
    }

    /// `ReportWrite::resolve_object`, for a hold.
    pub(crate) fn resolve(&self, pieces: &[usize]) -> std::io::Result<()> {
        self.resolve_object(pieces)
    }

    /// `~CBlastFormat`'s reset of the stream's exceptions, before the writes of a format that
    /// NCBI makes after it (SAM's `~CSAM_Formatter`, a member of `CBlastFormat`).
    pub fn stop_throwing(&self) {
        lock(&self.0).throws = false;
    }

    /// `ReportStream::finish`.
    pub fn finish(&self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        lock(&self.0).finish(result)
    }
}

/// A writer of a report that declares NCBI's regions and flushes to the command line's report
/// stream behind it (`ReportStream`, through `NcbiOstream`); every other writer ignores them.
pub trait ReportWrite: std::io::Write {
    /// The insertions from here on are made in NCBI's region `site`.
    fn site(&mut self, _site: Site) -> std::io::Result<()> {
        Ok(())
    }

    /// A flush of NCBI's (`flush`, `endl`): NCBI's buffer writes what it holds.
    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }

    /// After `is_undecided_report_write_failure`: NCBI's pieces of the object that LOSAT inserted
    /// whole (`Site::Serial { whole: true }`), which decide NCBI's end (an error: the run stops).
    fn resolve_object(&mut self, _pieces: &[usize]) -> std::io::Result<()> {
        Ok(())
    }

    /// An insertion of the report that sets the stream's `badbit` itself (a null `const char*`: the
    /// epilog's matrix name, `report::pairwise::write_matrix_line`), which throws on NCBI's stream
    /// (its `exceptions` mask holds `NcbiBadbit`, blast_format.cpp:118-119): the bytes the stream
    /// holds are written, the run stops, and the application reports "BLAST failed to write
    /// output", exit status 6 (blast_app_util.hpp:252-255). A writer that is no stream of NCBI's
    /// gives the error of a failed write.
    fn badbit(&mut self) -> std::io::Result<()> {
        Err(write_failure())
    }

    /// After `is_undecided_report_write_failure` of an object written once the formatter ended
    /// (the error archive): NCBI's write fault at `position` (`ReportStream`'s 802.5 post).
    fn write_fault(&mut self, _position: usize) -> std::io::Result<()> {
        Ok(())
    }

    /// The insertions `report[from..at]` for each cut `(at, _)` of `cuts`, then `report[..to]`
    /// from the last cut: one region's run of NCBI's insertions (every cut a `Cut::Insertion`,
    /// in increasing order, none past `to`), each one write of the stream (`CutCursor::write`).
    fn insert_cut_run(
        &mut self,
        report: &[u8],
        from: usize,
        cuts: &[(usize, Cut)],
        to: usize,
    ) -> std::io::Result<()> {
        let mut at = from;
        for end in cuts.iter().map(|&(end, _)| end).chain([to]) {
            let end = end.max(at);
            if end > at {
                self.write_all(&report[at..end])?;
            }
            at = end;
        }
        Ok(())
    }
}

/// The command line's report stream itself (after the formatter: the error archive).
impl ReportWrite for ReportStream {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        lock(&self.0).set_site(site)
    }

    fn insert_cut_run(
        &mut self,
        report: &[u8],
        from: usize,
        cuts: &[(usize, Cut)],
        to: usize,
    ) -> std::io::Result<()> {
        lock(&self.0).insert_cut_run(report, from, cuts, to)
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        lock(&self.0).flush()
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        lock(&self.0).resolve_object(pieces)
    }

    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        lock(&self.0).write_fault(position)
    }

    fn badbit(&mut self) -> std::io::Result<()> {
        lock(&self.0).badbit()
    }
}

impl<W: ReportWrite + ?Sized> ReportWrite for &mut W {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        (**self).site(site)
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        (**self).ncbi_flush()
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        (**self).resolve_object(pieces)
    }

    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        (**self).write_fault(position)
    }

    fn badbit(&mut self) -> std::io::Result<()> {
        (**self).badbit()
    }

    fn insert_cut_run(
        &mut self,
        report: &[u8],
        from: usize,
        cuts: &[(usize, Cut)],
        to: usize,
    ) -> std::io::Result<()> {
        (**self).insert_cut_run(report, from, cuts, to)
    }
}

/// LOSAT's buffer before the stream: it writes what it holds first.
impl<W: ReportWrite> ReportWrite for std::io::BufWriter<W> {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        std::io::Write::flush(self)?;
        self.get_mut().site(site)
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        std::io::Write::flush(self)?;
        self.get_mut().ncbi_flush()
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        self.get_mut().resolve_object(pieces)
    }

    fn badbit(&mut self) -> std::io::Result<()> {
        std::io::Write::flush(self)?;
        self.get_mut().badbit()
    }
}

/// A buffer of the report's bytes (a piece formatted before it is written).
impl ReportWrite for Vec<u8> {}

/// What NCBI's formatter does at an offset of a part of the report that LOSAT formats into a
/// buffer before it writes it (`Cuts`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cut {
    /// One insertion of NCBI's ends here and the next starts.
    Insertion,
    /// The insertions from here on are made in NCBI's region `site` (`ReportWrite::site`); the
    /// insertion before ends here.
    Site(Site),
    /// NCBI flushes the stream here (`endl`, `ReportWrite::ncbi_flush`); the insertion before ends
    /// here.
    Flush,
}

/// The places in a part of the report formatted into a buffer where NCBI's insertions end and
/// where NCBI's formatter enters a region or flushes (`Cut`), at offsets of the buffer in
/// increasing order: the buffer is written as NCBI inserts it (`CutCursor::write`), each insertion
/// one write of the stream (`insert`), so that the stream's buffer (`ReportStream`) writes and
/// fails where NCBI's does. A buffer without cuts is one insertion.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Cuts(Vec<(usize, Cut)>);

impl Cuts {
    /// `cut` at `at`, after the cuts so far.
    pub fn push(&mut self, at: usize, cut: Cut) {
        debug_assert!(self.0.last().is_none_or(|&(last, _)| last <= at));
        self.0.push((at, cut));
    }

    /// Room for `additional` more cuts (one allocation for a document's cuts).
    pub fn reserve(&mut self, additional: usize) {
        self.0.reserve(additional);
    }

    /// The cuts of `other`, after these.
    pub fn extend(&mut self, other: Cuts) {
        if self.0.is_empty() {
            self.0 = other.0;
            return;
        }
        debug_assert!(other
            .0
            .first()
            .is_none_or(|&(first, _)| { self.0.last().is_none_or(|&(last, _)| last <= first) }));
        self.0.extend(other.0);
    }

    /// A cursor that writes the buffer from its start.
    pub fn cursor(&self) -> CutCursor<'_> {
        CutCursor {
            cuts: &self.0,
            next: 0,
            region: (0, Site::Report),
            skipping: false,
            caught: false,
        }
    }
}

/// A part of the report formatted into a buffer (`bytes`) away from the stream (on a pool's
/// thread), with NCBI's regions and flushes made there (`cuts`), to be written with them
/// (`CutCursor::write`).
#[derive(Default)]
pub struct CutBuffer {
    /// The part's bytes.
    pub bytes: Vec<u8>,
    /// The part's regions and flushes.
    pub cuts: Cuts,
}

impl std::io::Write for CutBuffer {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.bytes.extend_from_slice(buf);
        Ok(buf.len())
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        self.bytes.extend_from_slice(buf);
        Ok(())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl ReportWrite for CutBuffer {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        self.cuts.push(self.bytes.len(), Cut::Site(site));
        Ok(())
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        self.cuts.push(self.bytes.len(), Cut::Flush);
        Ok(())
    }
}

/// NCBI's insertions of `report`, one query's report of outfmt 16 or 15 (`site`,
/// `Site::Xml2Report`), into the report's stream: the writes of its object stream's buffer, taken
/// back from the report's bytes (`out_buffer::xml_pieces`, `json_pieces`).
pub(crate) fn report_pieces(site: Site, report: &[u8]) -> Vec<usize> {
    match site {
        Site::Xml2Report { json: true } => crate::report::serial::out_buffer::json_pieces(report),
        _ => crate::report::serial::out_buffer::xml_pieces(report, |tag| {
            crate::report::blast_xml2::is_string_element(tag)
        }),
    }
}

/// The writing of a buffer with its `Cuts`, part by part (`write`).
pub struct CutCursor<'c> {
    cuts: &'c [(usize, Cut)],
    next: usize,
    /// Where the region being written starts, and its site.
    region: (usize, Site),
    /// The rest of the region was done where NCBI's pieces decided its end (`failed`), or where
    /// NCBI's formatter caught its failure.
    skipping: bool,
    /// A failure that NCBI's formatter catches ended a region in the part being written: the
    /// part's write returns it, after what follows the region (`write`).
    caught: bool,
}

impl CutCursor<'_> {
    /// Writes `report[from..to]` to `writer`, the part after the parts written before: an
    /// insertion up to each cut in it, and each cut's region or flush (those at `to` too, which
    /// come before what follows the part).
    pub fn write(
        &mut self,
        writer: &mut dyn ReportWrite,
        report: &[u8],
        from: usize,
        to: usize,
    ) -> std::io::Result<()> {
        self.write_until(writer, report, from, to, usize::MAX)
    }

    /// `write` with the cuts before the `cuts_end`th alone: those after it, at `to`, come after
    /// what is written between (a text on standard error that the formatter made after them,
    /// `blastx::blast_xml::HeldOutput`).
    pub fn write_until(
        &mut self,
        writer: &mut dyn ReportWrite,
        report: &[u8],
        from: usize,
        to: usize,
        cuts_end: usize,
    ) -> std::io::Result<()> {
        let cuts = &self.cuts[..cuts_end.min(self.cuts.len())];
        let mut at = from;
        loop {
            // The run of insertions up to the next region or flush (or `to`), in one call.
            let run = self.next;
            let mut next = run;
            while cuts
                .get(next)
                .is_some_and(|&(cut_at, cut)| cut_at <= to && cut == Cut::Insertion)
            {
                next += 1;
            }
            let other = cuts.get(next).copied().filter(|&(cut_at, _)| cut_at <= to);
            let stop = other.map_or(to, |(cut_at, _)| cut_at).max(at);
            if !self.skipping && (stop > at || next > run) {
                if let Err(error) = writer.insert_cut_run(report, at, &cuts[run..next], stop) {
                    self.failed(writer, report, error)?;
                }
            }
            self.next = next;
            at = stop;
            let Some((cut_at, cut)) = other else {
                // The catch returns to the formatter after the part (BLAST XML's next query, or
                // what follows a held XML2 report: its next insertion into the bad stream throws).
                if std::mem::take(&mut self.caught) {
                    return Err(caught_failure());
                }
                return Ok(());
            };
            self.next += 1;
            if self.skipping && !matches!(cut, Cut::Site(_)) {
                continue;
            }
            self.skipping = false;
            match cut {
                Cut::Insertion => {}
                Cut::Site(site) => {
                    if let Err(error) = writer.site(site) {
                        // The region's writes, made at its end, failed where NCBI's formatter
                        // catches it: that region is done, and the new one has started.
                        self.region_failed(writer, report, error)?;
                        self.skipping = false;
                    }
                    self.region = (cut_at, site);
                }
                Cut::Flush => {
                    if let Err(error) = writer.ncbi_flush() {
                        self.failed(writer, report, error)?;
                    }
                }
            }
        }
    }

    /// Whether the rest of the region being written is not written (`failed`).
    pub fn skipping(&self) -> bool {
        self.skipping
    }

    /// A failure of the region being written outside its bytes (the tie of a text on standard
    /// error that NCBI's formatter catches, `is_caught_report_write_failure`): the rest of the
    /// region is skipped, and the formatter's next insertions are in `Site::Report` (`failed`).
    pub fn caught_outside(&mut self, writer: &mut dyn ReportWrite) -> std::io::Result<()> {
        self.skipping = true;
        writer.site(Site::Report)
    }

    /// A failed insertion or flush of the region being written (`region_failed`). Where NCBI's
    /// formatter catches it, the catch returns to the formatter, whose next insertions are in
    /// `Site::Report` until it declares another region.
    #[cold]
    fn failed(
        &mut self,
        writer: &mut dyn ReportWrite,
        report: &[u8],
        error: std::io::Error,
    ) -> std::io::Result<()> {
        if self.region_failed(writer, report, error)? {
            writer.site(Site::Report)?;
        }
        Ok(())
    }

    /// A failed write of the region being written. A report of outfmt 16 or 15 is inserted
    /// whole (`Site::Xml2Report`): where its write fails at the device, NCBI's pieces of the
    /// report (`out_buffer::xml_pieces`, `json_pieces`, from its bytes) decide NCBI's end
    /// (`ReportWrite::resolve_object`); where only its `Close` fails, the rest of the region
    /// (`Close`'s flush) is done.
    /// A failure that NCBI's formatter catches (`is_caught_report_write_failure`) ends the region:
    /// the rest of it is skipped, and what follows it is written (`caught`); returns whether it
    /// was caught.
    #[cold]
    fn region_failed(
        &mut self,
        writer: &mut dyn ReportWrite,
        report: &[u8],
        error: std::io::Error,
    ) -> std::io::Result<bool> {
        let resolved = if is_undecided_report_write_failure(&error) {
            let (start, site) = self.region;
            let end = self.cuts[self.next..]
                .iter()
                .find(|(_, cut)| matches!(cut, Cut::Site(_)))
                .map_or(report.len(), |&(at, _)| at);
            let pieces = report_pieces(site, &report[start.min(end)..end]);
            writer.resolve_object(&pieces)
        } else {
            Err(error)
        };
        let caught = match resolved {
            Err(error) if is_caught_report_write_failure(&error) => true,
            resolved => {
                resolved?;
                false
            }
        };
        self.caught |= caught;
        self.skipping = true;
        Ok(caught)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1911-1931
/// NCBI's `DisplaySeqalign` of a query's alignments: from its first alignment (none shown: no
/// call that writes) they are written in its region (`Site::Alignments`), and at their end the
/// destructor of its `out2` flushes the stream (`CObjectOStream::Close`, `objostr.cpp:338-354`),
/// before the query's footer.
#[derive(Default)]
pub(crate) struct DisplaySeqalign {
    open: bool,
}

impl DisplaySeqalign {
    /// Before an alignment: the region starts at the first.
    pub(crate) fn enter<W: ReportWrite + ?Sized>(&mut self, writer: &mut W) -> std::io::Result<()> {
        if !self.open {
            self.open = true;
            writer.site(Site::Alignments)?;
        }
        Ok(())
    }

    /// After the query's last alignment: `out2`'s flush, then the report's region.
    pub(crate) fn leave<W: ReportWrite + ?Sized>(&mut self, writer: &mut W) -> std::io::Result<()> {
        if self.open {
            self.open = false;
            writer.ncbi_flush()?;
            writer.site(Site::Report)?;
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:456
/// An insertion of `NcbiEndl` (`endl`): a line feed, then NCBI's flush (`ReportWrite::ncbi_flush`),
/// as at the start of each query's footer and twice at the start of a pairwise report's epilog.
pub fn ncbi_endl<W: ReportWrite + ?Sized>(out: &mut W) -> std::io::Result<()> {
    out.write_all(b"\n")?;
    out.ncbi_flush()
}

/// A stream insertion (`operator<<`) in NCBI's report: one insertion (a string or a number) of
/// `bytes`, one write call, which the command line's report stream takes as one (`ReportStream`);
/// other writers just take the bytes. A string that LOSAT streams is `insert_streamed`.
pub fn insert(out: &mut (impl std::io::Write + ?Sized), bytes: &[u8]) -> std::io::Result<()> {
    out.write_all(bytes)
}

/// One insertion of NCBI's of `len` bytes that LOSAT streams instead of holding (the organism
/// report, whose runs of about 4 GiB of spaces `tax_report` keeps as counts): `fill` gives the bytes to
/// the `Insertion`, which writes them in pieces that NCBI's buffers (`NcbiBuffer`) write as they
/// write the whole. Each piece is longer than a buffer (an output file's buffer writes it with
/// what it holds and keeps nothing; standard output's keeps `(held + n) % block`, which adds up
/// over the pieces as over the whole, `cli::NcbiBuffer::sputn`), and each but the last ends in `NCBI_BUFSIZ` bytes without a newline (so a
/// line-buffered `cout` keeps its whole remainder there, as within the whole); the last is the
/// rest. A string without such a place is written in one piece.
pub fn insert_streamed<W: std::io::Write + ?Sized>(
    out: &mut W,
    len: u64,
    fill: impl FnOnce(&mut Insertion<'_, W>) -> std::io::Result<()>,
) -> std::io::Result<()> {
    let mut insertion = Insertion {
        out,
        piece: Vec::new(),
        left: len,
        since_newline: 0,
    };
    fill(&mut insertion)?;
    debug_assert_eq!(insertion.left, insertion.piece.len() as u64);
    insertion.out.write_all(&insertion.piece)
}

/// The bytes of one insertion that LOSAT streams (`insert_streamed`).
pub struct Insertion<'w, W: std::io::Write + ?Sized> {
    out: &'w mut W,
    /// The piece being made.
    piece: Vec<u8>,
    /// The bytes of the insertion not written yet, `piece`'s among them.
    left: u64,
    /// How many bytes follow the last newline.
    since_newline: u64,
}

impl<W: std::io::Write + ?Sized> Insertion<'_, W> {
    /// The length from which a piece is written where it may end.
    const PIECE: usize = 1 << 16;

    /// The next bytes of the insertion.
    pub fn bytes(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        self.piece.extend_from_slice(bytes);
        self.since_newline = match bytes.iter().rposition(|&byte| byte == b'\n') {
            Some(at) => (bytes.len() - at - 1) as u64,
            None => self.since_newline + bytes.len() as u64,
        };
        self.end_piece()
    }

    /// The next `count` bytes of the insertion, spaces.
    pub fn spaces(&mut self, mut count: u64) -> std::io::Result<()> {
        while count > 0 {
            let n = count.min(Self::PIECE as u64);
            self.piece.resize(self.piece.len() + n as usize, b' ');
            self.since_newline += n;
            count -= n;
            self.end_piece()?;
        }
        Ok(())
    }

    /// Writes the piece where it may end (`insert_streamed`) once it is `PIECE` bytes long. A
    /// piece ends only if more than `NCBI_BUFSIZ` bytes follow it, which is at least a block (a
    /// block never exceeds `NCBI_BUFSIZ`), so the last piece always takes the buffer's branch of
    /// an insertion longer than the room, as the whole insertion does; a final piece that exactly
    /// filled a block would differ from the whole.
    fn end_piece(&mut self) -> std::io::Result<()> {
        let n = self.piece.len() as u64;
        if n >= Self::PIECE as u64
            && self.since_newline >= NCBI_BUFSIZ as u64
            && self.left - n > NCBI_BUFSIZ as u64
        {
            self.out.write_all(&self.piece)?;
            self.left -= n;
            self.piece.clear();
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:68-93
/// The stream of a report as NCBI's formatter writes to it (`CBlastFormat::m_Outfile`), with
/// nothing of a program's before it: each write is one insertion (`operator<<`) of NCBI's
/// (`insert`, `insert_streamed`), passed on as one write. A writer becomes one only through the
/// constructors below, one per destination, each saying what comes after it: they are the only
/// places where a buffer or hold could come between the report's insertions and NCBI's buffer
/// (`ReportStream`, which counts each write as one insertion). A format's sink (`sink`: the
/// command line's `ReportStream`, `OutputSink::Report`, the same stream of the library's
/// `run_to_file`, or the library's and LOSAT Web's writers), blastx's command line (`report`), and
/// the hold of a library call (`held`). Only the command line's `-out` and standard output streams
/// and the `-out` stream of `run_to_file` take NCBI's regions and flushes (`ReportWrite`).
pub struct NcbiOstream<'a> {
    out: Box<dyn std::io::Write + 'a>,
    report: Option<ReportHandle>,
    /// The report stream that `out` is (`report`'s), written without `out`.
    direct: Option<ReportHandle>,
    /// A hold before the report stream (`before_report`), which takes the regions and flushes in
    /// order with its bytes.
    hold: Option<Box<dyn ReportWrite + 'a>>,
}

impl<'a> NcbiOstream<'a> {
    fn new(out: Box<dyn std::io::Write + 'a>) -> Self {
        NcbiOstream {
            out,
            report: None,
            direct: None,
            hold: None,
        }
    }

    /// A format's sink opened as its report's stream (`OutputSink::open`): the command line's
    /// `ReportStream` (`OutputSink::Report`) and the same stream of the library's `run_to_file`
    /// (`OutputSink::Out`, `ReportStream::of_library_file`); or the caller's writer as it is,
    /// which a library caller or LOSAT Web gives (no stream of NCBI's to model); or, for a
    /// library call's `Stdout` and `File` sinks, LOSAT's buffered standard output or file (the
    /// file kept for the run, `api::local_blast::KeptFile`; no stream of NCBI's either: the
    /// report's bytes are NCBI's, the times at which they reach the device are not).
    pub(crate) fn sink(
        sink: &'a mut crate::api::local_blast::OutputSink<'_>,
    ) -> std::io::Result<Self> {
        let report = sink.report_handle();
        let out = sink.open()?;
        Ok(NcbiOstream {
            out,
            direct: report.clone(),
            report,
            hold: None,
        })
    }

    /// The command line's report stream `stream` of blastx and rpstblastn (`native::run_inner`),
    /// every format's.
    pub(crate) fn report(stream: ReportStream) -> NcbiOstream<'static> {
        let report = stream.handle();
        NcbiOstream {
            out: Box::new(stream),
            direct: Some(report.clone()),
            report: Some(report),
            hold: None,
        }
    }

    /// `hold`, a hold of blastx's command line before its report stream `report`
    /// (`blast_xml::HeldOutput`): it takes the bytes, and NCBI's regions and flushes in order with
    /// them (it passes them to `report`).
    pub(crate) fn before_report<H: ReportWrite + Clone + 'a>(
        hold: H,
        report: ReportHandle,
    ) -> Self {
        NcbiOstream {
            out: Box::new(hold.clone()),
            report: Some(report),
            direct: None,
            hold: Some(Box::new(hold)),
        }
    }

    /// The command line's report stream behind this one (`report`, `sink`).
    pub(crate) fn report_handle(&self) -> Option<&ReportHandle> {
        self.report.as_ref()
    }

    /// A format's stream (`sink`'s) behind the hold of a library call (`local::report_batches`,
    /// `blast_xml::HeldFormats`), which passes the writes and flushes on in their order (at the
    /// release, while it holds). The caller's writer has no stream of NCBI's (`sink`), so no buffer
    /// of NCBI's comes after the hold.
    pub(crate) fn held(
        stream: crate::algorithm::blastx::HeldFormatWriter<NcbiOstream<'a>>,
    ) -> Self {
        Self::new(Box::new(stream))
    }
}

impl std::io::Write for NcbiOstream<'_> {
    /// One insertion, written with one write of the stream.
    #[inline]
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.write_all(buf)?;
        Ok(buf.len())
    }

    #[inline]
    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        match &self.direct {
            Some(report) => report.insert(buf),
            None => self.out.write_all(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match &self.direct {
            // LOSAT's flushes of the command line's stream hold (`ReportStream`'s).
            Some(_) => Ok(()),
            None => self.out.flush(),
        }
    }
}

impl ReportWrite for NcbiOstream<'_> {
    fn site(&mut self, site: Site) -> std::io::Result<()> {
        if let Some(hold) = &mut self.hold {
            return hold.site(site);
        }
        match &self.report {
            Some(report) => report.site(site),
            None => Ok(()),
        }
    }

    fn insert_cut_run(
        &mut self,
        report: &[u8],
        from: usize,
        cuts: &[(usize, Cut)],
        to: usize,
    ) -> std::io::Result<()> {
        if let Some(hold) = &mut self.hold {
            return hold.insert_cut_run(report, from, cuts, to);
        }
        match &self.direct {
            Some(stream) => lock(&stream.0).insert_cut_run(report, from, cuts, to),
            None => self.out.write_all(&report[from..to.max(from)]),
        }
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        if let Some(hold) = &mut self.hold {
            return hold.ncbi_flush();
        }
        match &self.report {
            Some(report) => report.flush(),
            None => Ok(()),
        }
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        match &self.report {
            Some(report) => report.resolve_object(pieces),
            None => Ok(()),
        }
    }

    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        match &self.report {
            Some(report) => report.write_fault(position),
            None => Ok(()),
        }
    }

    fn badbit(&mut self) -> std::io::Result<()> {
        match (&self.report, &mut self.hold) {
            (_, Some(hold)) => hold.badbit(),
            (Some(report), None) => report.badbit(),
            (None, None) => Err(write_failure()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:768-778
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:893-896 (`COStreamBuffer::FlushBuffer`)
// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:339-354, :385-395
/// Writes `document`, the report of a side file of outfmt 13 or 14, to `file` as NCBI's
/// `CNcbiOfstream` (no exception mask; its file buffer, `NcbiStdio::File`) gets it from the buffer
/// of its object stream: `pieces` are that buffer's writes (`out_buffer::xml_pieces`,
/// `json_pieces`; `back_limit` its `m_BackLimit`), each one insertion of the stream.
///
/// A piece whose write of the file fails is the buffer's "write fault" (`CIOException::eWrite`):
/// the stream's `Close` posts 802.5 at the buffer's position then (the bytes it had written and
/// those it held, the piece's and the back limit's) with an empty stack, and the exception ends
/// the run at `CATCH_ALL`, which swallows that code: exit status 0, no later report and no master
/// file (an error of exit status 0, its text posted here). Where only `Close`'s flush fails, its
/// chain of posts ends in the destructor's 802.4 and the run goes on. Without the command line's
/// run (a library call writes its side files without a stream of NCBI's) a failed write is the
/// old error, `blast_xml2::framing::side_file_write_error`.
pub fn write_side_file(
    mut file: std::fs::File,
    document: &[u8],
    pieces: impl FnOnce() -> Vec<usize>,
    back_limit: usize,
) -> anyhow::Result<()> {
    use std::io::Write;
    // NCBI's writes of the report before the side file are made first.
    sync_report()?;
    let Some(program) = *lock(&RUN_PROGRAM) else {
        return file
            .write_all(document)
            .and_then(|()| file.flush())
            .map_err(|_| crate::report::blast_xml2::framing::side_file_write_error());
    };
    // The document is written whole; where the file fails, NCBI's pieces decide where NCBI's
    // writes failed (`side_file_failure`).
    let mut landed = 0u64;
    if write_device(&mut file, &mut landed, document)
        .and_then(|()| file.flush())
        .is_ok()
    {
        return Ok(());
    }
    side_file_failure(program, landed, &pieces(), back_limit)
}

/// NCBI's end where a side file's device took `landed` bytes of its document, whose object
/// stream's buffer inserts `pieces` into the file's buffer (`NcbiStdio::File`, `write_side_file`).
#[cold]
#[inline(never)]
fn side_file_failure(
    program: &str,
    landed: u64,
    pieces: &[usize],
    back_limit: usize,
) -> anyhow::Result<()> {
    use std::io::Write;
    // A write of the file's buffer inside the object that passes the bytes the file took.
    let buffer = NcbiBuffer::new(NcbiStdio::File);
    if let Some((k, end)) = failed_piece(buffer, 0, pieces, landed) {
        // `EndOfWrite`'s write leaves nothing in the object stream's buffer.
        let leave = if k + 1 == pieces.len() { 0 } else { back_limit };
        let post = format!(
            "Error: (802.5) [{program}] CObjectOStream: error at byte {}: stack is empty: \
             write fault\n",
            end as usize + leave
        );
        return Err(crate::cli::NativeError {
            exit: 0,
            message: post.into_bytes(),
        }
        .into());
    }
    let end: usize = pieces.iter().sum();
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:899-907 (`COStreamBuffer::Flush`)
    // NCBI's chain, written in LOSAT's words (`crate::wording::losat_bytes`: without the
    // exceptions' classes and codes, `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12) by the
    // run's standard error (`CliStderr`, tied to the report).
    let position = end;
    let chain = format!(
        "Error: (802.5) [{program}] CObjectOStream: error at byte {position}: stack is empty: \
         cannot close output stream\n\
         Error: ({FLUSH_FAILED}) [{program}] COStreamBuffer::Flush: failed\n\
         Error: ({SERIAL_IO_ERROR}) [{program}] byte {position}: cannot close output stream\n\
         Error: (802.4) [{program}] Cannot close serializing output stream \
         ({SERIAL_IO_ERROR}) byte {position}: cannot close output stream\n"
    );
    let mut stderr = CliStderr::process();
    stderr.write_all(&crate::wording::losat_bytes(chain.as_bytes()))?;
    stderr.flush()?;
    Ok(())
}

/// `GetType()::GetErrCodeString()` of the failed flush of a side file's object stream.
const FLUSH_FAILED: &str = "CIOException::eFlush";

/// `GetType()::GetErrCodeString()` of the serial stream's error of a failed close.
const SERIAL_IO_ERROR: &str = "CSerialException::eIoError";

/// `cerr.tie()`: before a write to standard error, the run's report on standard output writes
/// what it holds (`ReportStream::standard_output`); a report in `-out` makes NCBI's writes so far
/// (`Core::sync`).
pub fn tie_standard_output() -> std::io::Result<()> {
    match tied_report() {
        Some(report) => {
            let mut report = lock(&report);
            if report.standard_output {
                report.tie()
            } else {
                // The `-out` file is not tied: NCBI's writes before the text are made, where one
                // that failed ended NCBI's run before it.
                report.sync()
            }
        }
        None => Ok(()),
    }
}

/// The run's report stream (`TIED_REPORT`) while it lives.
fn tied_report() -> Option<std::sync::Arc<std::sync::Mutex<Core>>> {
    lock(&TIED_REPORT)
        .as_ref()
        .and_then(std::sync::Weak::upgrade)
}

/// Before a text that NCBI writes outside `cerr`, where no tie flushes `cout`: with C's `stderr`
/// (`c_stderr_write`), through the handler of the `-logfile` file (`cerr_post`) or to a side file
/// (`write_side_file`). NCBI's writes of the report so far are made (`Core::sync`); what NCBI's
/// buffer holds stays there.
pub fn sync_report() -> std::io::Result<()> {
    match tied_report() {
        Some(report) => lock(&report).sync(),
        None => Ok(()),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3410-3412
/// `text` on the command line's NCBI C `stderr`, which NCBI's formatter writes with `fprintf`
/// (`InitConfig`'s `REGISTRY:` line of `GETURL_DEBUG`), apart from its `cerr` (`CliStderr`) and
/// from the run's diagnostics: no tie flushes `cout`, NCBI's writes of the report so far are made
/// first (`sync_report`), so the line comes before the bytes that NCBI's buffer holds where both
/// streams go to one pipe (line RG T4a2 fix round 2, review c C4: after the prolog's last
/// `NcbiEndl`; review e S1: after a first query's footer `NcbiEndl` where that query has no hits,
/// probes `T/t4/review_t4a2_e/out/nh1.*`, `nhx.*`). The text goes to the terminal where the run's
/// messages follow `-logfile` (`ncbi_environment::terminal_stderr`: C's `stderr` stays the
/// terminal), and a failed write is dropped, also on a bad `cerr` (NCBI reads nothing that
/// `fprintf` returns). The one writer of that stream, for the command line's process only
/// (`query_warnings::QueryWarnings::html_stderr_line`): a library call has no C `stderr` of its
/// own, its line goes to the call's diagnostics.
pub(crate) fn c_stderr_write(text: &[u8]) -> std::io::Result<()> {
    use std::io::Write;
    if text.is_empty() {
        return Ok(());
    }
    sync_report()?;
    let _ = match crate::blastinput::ncbi_environment::terminal_stderr() {
        Some(mut terminal) => terminal.write_all(text).and_then(|()| terminal.flush()),
        None => {
            let mut stderr = std::io::stderr();
            stderr.write_all(text).and_then(|()| stderr.flush())
        }
    };
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:6695-6696
/// NCBI's `cerr` (unbuffered, tied to `cout`) as NCBI BLAST+ 2.17.0's Linux build treats a write
/// that fails (observed: standard error on `/dev/full`, in a file at its size limit, or a pipe
/// whose reader closed it; line RG item T6a, review 2 D2): the stream becomes bad and the run goes
/// on as if the text were written; a bad `cerr` takes no later text, and a later post no longer
/// makes standard output write what it holds first (the tie, `tie_standard_output`, works on a
/// good `cerr` only). So a full standard error changes neither the report nor the exit status. A
/// pipe whose reader closed it fails the same way (`EPIPE`: the Rust runtime ignores SIGPIPE),
/// where NCBI is ended by SIGPIPE (exit status 141, no text), which LOSAT cannot be without a
/// native call: the form of NCBI's failed write of standard error (approved exception 5 of
/// `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, the Owner's answer A to line RG's T6-Q1, 2026-11-17).
/// Where the standard error is the file of `-logfile` (tblastx's, rpsblast's and rpstblastn's
/// run, `DiagTarget::into_standard_error`), a failed write is dropped and the next is made, as
/// NCBI's handler of that file writes each post and ignores what the write returns
/// (`CFileHandleDiagHandler::Post`; NCBI BLAST+ 2.17.0 tblastx `-logfile /dev/full` runs and
/// writes its report; line WQA's D11 round 2, 2026-11-20).
pub(crate) struct Cerr(std::sync::atomic::AtomicBool);

/// The command line's `cerr`: the process's standard error.
static PROCESS_CERR: Cerr = Cerr::new();

impl Cerr {
    const fn new() -> Self {
        Cerr(std::sync::atomic::AtomicBool::new(false))
    }

    /// Whether a write failed (`badbit`).
    fn failed(&self) -> bool {
        self.0.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Writes `text` to `device`; a failed write makes the stream bad.
    fn write(&self, device: &mut dyn std::io::Write, text: &[u8]) {
        if device
            .write_all(text)
            .and_then(|()| device.flush())
            .is_err()
        {
            self.0.store(true, std::sync::atomic::Ordering::Relaxed);
        }
    }

    /// Writes `text` to the process's standard error (`write`), or to the file of `-logfile`,
    /// whose failed write is dropped.
    fn write_process(&self, text: &[u8]) {
        if crate::blastinput::ncbi_environment::terminal_stderr().is_some() {
            use std::io::Write;
            let mut stderr = std::io::stderr();
            let _ = stderr.write_all(text).and_then(|()| stderr.flush());
            return;
        }
        self.write(&mut std::io::stderr(), text);
    }
}

/// `text` on the command line's standard error as NCBI's `cerr` takes it (`Cerr`): after the tie
/// (`tie_standard_output`, whose error is the report stream's: `std::terminate` while the
/// formatter lives), nothing on a bad `cerr`. The application's posts outside the report (the
/// error of `CATCH_ALL`, `cli::exit_on_native_error`) are written so too. Where the run's
/// messages follow `-logfile` (`ncbi_environment::terminal_stderr`), NCBI's handler of that file
/// writes them, not `cerr`, so no tie flushes `cout`: NCBI's writes so far are made alone
/// (`sync_report`; NCBI BLAST+ 2.17.0, line RG T4a2 review e S1, `tblastx -html -logfile` with
/// both streams in one pipe, probes `T/t4/review_t4a2_e/out/txlog.rep.*`: `GETURL_DEBUG`'s
/// `REGISTRY:` line after the version line's flush, 166 bytes, where a flush at each post put it
/// after the first description table's head, 1,037).
pub fn cerr_post(text: &[u8]) -> std::io::Result<()> {
    if text.is_empty() || PROCESS_CERR.failed() {
        return Ok(());
    }
    if crate::blastinput::ncbi_environment::terminal_stderr().is_some() {
        sync_report()?;
    } else {
        tie_standard_output()?;
    }
    PROCESS_CERR.write_process(text);
    Ok(())
}

/// The command line's standard error (`cerr_post`), tied to the report's standard output
/// (`tie_standard_output`): what standard output writes before a text is written before it where
/// the text is made. A whole message is written with one `write_all` of the stream, which
/// `std::io::Stderr` writes under one lock of standard error, so that no other thread's message
/// comes inside it. An empty write is no write of NCBI's to `cerr` (a query without warnings,
/// whose `ERR_POST` NCBI does not make): it writes nothing of the report. A failed write of
/// standard error is no error of the run (`Cerr`).
pub struct CliStderr;

impl CliStderr {
    /// The process's standard error.
    pub fn process() -> Self {
        CliStderr
    }
}

impl std::io::Write for CliStderr {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        cerr_post(buf)?;
        Ok(buf.len())
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        cerr_post(buf)
    }

    fn write_fmt(&mut self, args: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        self.write_all(args.to_string().as_bytes())
    }

    /// Standard error holds nothing.
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[cfg(all(test, unix))]
#[path = "ostream_tests.rs"]
mod tests;
