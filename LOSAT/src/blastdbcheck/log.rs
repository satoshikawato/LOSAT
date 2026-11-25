//! blastdbcheck's report (`CBlastDbCheckLog`): the stream of `-logfile` or the standard output,
//! written at the verbosity of `-verbosity`, and what ends a run before NCBI's would.

use std::io::Write;

use crate::ncbi_app::files::Ostream;

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:95-103
/// The verbosity levels: a line of level `L` is written when `L <= -verbosity`.
pub(super) const BRIEF: i32 = 1;
pub(super) const SUMMARY: i32 = 2;
pub(super) const DETAILS: i32 = 3;
pub(super) const MINUTIAE: i32 = 4;

/// What ends a run before NCBI's run ends.
pub(super) enum Stop {
    /// A write found the report's pipe closed, where SIGPIPE ends NCBI's process
    /// (`files::Ostream::closed`): the run ends at once, with exit status 0 and no message.
    Closed,
    /// LOSAT's refusal: what NCBI does there is not what LOSAT gives (the reason).
    Refused(String),
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:258-287
// A line of a level above the verbosity goes to a string stream that is emptied at each use
// (`m_DevNull`): nothing. NCBI's lines end with `endl` (or `flush`), so each is written and
// flushed whole: the report's file of `-logfile` and the diagnostics' file of the same name take
// the writes in NCBI's order (`report::run_files`).
/// The report: its stream and the verbosity.
pub(super) struct Log<'a> {
    out: Ostream<'a>,
    level: i32,
}

impl<'a> Log<'a> {
    pub(super) fn new(out: &'a mut dyn Write, level: i32) -> Self {
        Log {
            out: Ostream::new(out),
            level,
        }
    }

    /// The verbosity (`GetLevel`).
    pub(super) fn level(&self) -> i32 {
        self.level
    }

    /// Whether a line of `level` is written.
    pub(super) fn shows(&self, level: i32) -> bool {
        level <= self.level
    }

    /// Writes `text` (whole lines with their line feeds) if a line of `level` is written.
    pub(super) fn write(&mut self, level: i32, text: &[u8]) -> Result<(), Stop> {
        if !self.shows(level) {
            return Ok(());
        }
        let _ = self.out.write_all(text);
        let _ = self.out.flush();
        if self.out.closed() {
            return Err(Stop::Closed);
        }
        Ok(())
    }

    /// Writes the line made of `parts` and a line feed, if a line of `level` is written.
    pub(super) fn line(&mut self, level: i32, parts: &[&[u8]]) -> Result<(), Stop> {
        if !self.shows(level) {
            return Ok(());
        }
        let mut text = parts.concat();
        text.push(b'\n');
        self.write(level, &text)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:983-984
    /// The text of an exception after a line that says it was caught: `<what()>` and the line feed
    /// of `endl` (a `CException`'s text ends with its own, so an empty line follows it), after
    /// `prefix`, at level Detailed. `what` is `None` where LOSAT does not give NCBI's text of that
    /// error: LOSAT's refusal, where the line is written.
    pub(super) fn what(
        &mut self,
        prefix: &[u8],
        what: Option<&[u8]>,
        error: &dyn Fn() -> String,
    ) -> Result<(), Stop> {
        if !self.shows(DETAILS) {
            return Ok(());
        }
        let Some(what) = what else {
            return Err(Stop::Refused(error()));
        };
        self.line(DETAILS, &[prefix, what])
    }
}
