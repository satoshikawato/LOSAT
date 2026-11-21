//! No NCBI counterpart: the command line's progress bar (`batch_progress_bar`) on a terminal
//! kept in memory, which stands for the standard error on a terminal.

use super::batch_progress_bar;
use indicatif::{ProgressDrawTarget, TermLike};
use std::io;
use std::sync::{Arc, Mutex};

/// The rows of a terminal of 80 columns and its cursor (row, column).
#[derive(Debug, Default)]
struct Screen {
    rows: Vec<Vec<char>>,
    row: usize,
    col: usize,
}

impl Screen {
    fn put(&mut self, ch: char) {
        while self.rows.len() <= self.row {
            self.rows.push(Vec::new());
        }
        let line = &mut self.rows[self.row];
        while line.len() <= self.col {
            line.push(' ');
        }
        line[self.col] = ch;
        self.col += 1;
    }
}

/// A terminal kept in memory; a clone writes to the same screen.
#[derive(Debug, Clone, Default)]
struct MemoryTerm(Arc<Mutex<Screen>>);

impl MemoryTerm {
    /// The rows that show something (trailing blanks dropped).
    fn shown(&self) -> Vec<String> {
        let screen = self.0.lock().unwrap();
        screen
            .rows
            .iter()
            .map(|row| row.iter().collect::<String>().trim_end().to_string())
            .filter(|row| !row.is_empty())
            .collect()
    }

    fn cursor(&self) -> (usize, usize) {
        let screen = self.0.lock().unwrap();
        (screen.row, screen.col)
    }
}

impl TermLike for MemoryTerm {
    fn width(&self) -> u16 {
        80
    }

    fn move_cursor_up(&self, n: usize) -> io::Result<()> {
        let mut screen = self.0.lock().unwrap();
        screen.row = screen.row.saturating_sub(n);
        Ok(())
    }

    fn move_cursor_down(&self, n: usize) -> io::Result<()> {
        self.0.lock().unwrap().row += n;
        Ok(())
    }

    fn move_cursor_right(&self, n: usize) -> io::Result<()> {
        self.0.lock().unwrap().col += n;
        Ok(())
    }

    fn move_cursor_left(&self, n: usize) -> io::Result<()> {
        let mut screen = self.0.lock().unwrap();
        screen.col = screen.col.saturating_sub(n);
        Ok(())
    }

    fn write_line(&self, s: &str) -> io::Result<()> {
        self.write_str(s)?;
        self.write_str("\n")
    }

    fn write_str(&self, s: &str) -> io::Result<()> {
        let mut screen = self.0.lock().unwrap();
        for ch in s.chars() {
            match ch {
                '\r' => screen.col = 0,
                '\n' => {
                    screen.row += 1;
                    screen.col = 0;
                }
                _ => screen.put(ch),
            }
        }
        Ok(())
    }

    fn clear_line(&self) -> io::Result<()> {
        let mut screen = self.0.lock().unwrap();
        let row = screen.row;
        if let Some(line) = screen.rows.get_mut(row) {
            line.clear();
        }
        screen.col = 0;
        Ok(())
    }

    fn flush(&self) -> io::Result<()> {
        Ok(())
    }
}

/// A bar drawn on `term` after three of five subjects.
fn drawn_bar(term: &MemoryTerm) -> indicatif::ProgressBar {
    let bar = batch_progress_bar(5, ProgressDrawTarget::term_like(Box::new(term.clone())));
    bar.inc(3);
    bar.tick();
    let shown = term.shown();
    assert_eq!(shown.len(), 1, "the bar is one line: {shown:?}");
    assert!(shown[0].contains("3/5"), "the bar shows 3/5: {shown:?}");
    bar
}

#[test]
fn the_end_of_a_batch_clears_the_bar() {
    let term = MemoryTerm::default();
    let bar = drawn_bar(&term);
    // What `search_query_batch` calls when the batch's search ends.
    bar.finish_using_style();
    assert!(bar.is_finished());
    assert_eq!(
        term.shown(),
        Vec::<String>::new(),
        "nothing of the bar stays"
    );
    assert_eq!(
        term.cursor(),
        (0, 0),
        "the next output starts where the bar started"
    );
    // The drop after the end writes nothing more.
    drop(bar);
    assert_eq!(term.shown(), Vec::<String>::new());
}

#[test]
fn a_search_that_stops_with_an_error_clears_the_bar() {
    let term = MemoryTerm::default();
    // A batch that returns an error drops its bar unfinished.
    drop(drawn_bar(&term));
    assert_eq!(
        term.shown(),
        Vec::<String>::new(),
        "nothing of the bar stays"
    );
    assert_eq!(term.cursor(), (0, 0));
}
