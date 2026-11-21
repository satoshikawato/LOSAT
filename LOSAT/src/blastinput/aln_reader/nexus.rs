//! NCBI's NEXUS scanner (`CAlnScannerNexus`): `#NEXUS`, then commands that end in `;`, each in a
//! block (`BEGIN name; ... END;`). The `DATA` and `CHARACTERS` blocks give the counts
//! (`DIMENSIONS`), the gap, missing and match characters (`FORMAT`) and the rows (`MATRIX`); the
//! `TAXA` block the count; the `NCBI` block the deflines (`SEQUIN`). Every other command and block
//! is read over. `CAlnReader` gives the scanner the sequence info of its own `SetClustal` call, so
//! the `DATATYPE` of the file changes nothing (the constructor's `m_UseNexusInfo` is not read by
//! `ReadAlignmentFile`).
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_nexus.cpp:61-960

use super::base::{check_id, process_defline, verify_rows, Scan, ScanFailure, Scanned, SeqInfo};
use super::lines::PeekAheadStream;
use crate::ncbi_str::{is_blank, string_to_int, truncate_spaces};
use crate::utils::format_guess::tokens;

/// The cause of the refusal of a command that has no text where NCBI's code reads its last line
/// (`list::back()` and `front()` of an empty list).
pub(crate) const EMPTY_COMMAND: &str =
    "NCBI BLAST+ 2.17.0 reads a NEXUS command that has no text through an empty list";

/// `SLineInfo`: a piece of a command and the line it came from.
#[derive(Clone)]
struct Line {
    data: Vec<u8>,
    num: i32,
}

type Args = Vec<Line>;

fn find_first_of(hay: &[u8], set: &[u8], from: usize) -> Option<usize> {
    hay.get(from..)?
        .iter()
        .position(|c| set.contains(c))
        .map(|at| at + from)
}

fn find_first_not_of(hay: &[u8], set: &[u8], from: usize) -> Option<usize> {
    hay.get(from..)?
        .iter()
        .position(|c| !set.contains(c))
        .map(|at| at + from)
}

/// `find_last_not_of(set, pos)`: the last position at or before `pos` whose byte is not in `set`.
fn find_last_not_of(hay: &[u8], set: &[u8], pos: usize) -> Option<usize> {
    if hay.is_empty() {
        return None;
    }
    let end = pos.min(hay.len() - 1);
    (0..=end).rev().find(|&i| !set.contains(&hay[i]))
}

/// `NStr::FindNoCase(str, pattern)`.
fn find_no_case(hay: &[u8], pattern: &[u8]) -> Option<usize> {
    let hay = hay.to_ascii_lowercase();
    let pattern = pattern.to_ascii_lowercase();
    if pattern.is_empty() {
        return Some(0);
    }
    hay.windows(pattern.len())
        .position(|w| w == pattern.as_slice())
}

/// `sFindCharOutsideComment`: the first `c` at or after `start` while no comment is open; the
/// brackets on the way change the open count.
fn find_char_outside_comment(c: u8, line: &[u8], open: &mut i32, start: usize) -> Option<usize> {
    for (index, &b) in line.iter().enumerate().skip(start) {
        if b == b'[' {
            *open += 1;
        } else if b == b']' {
            *open -= 1;
        } else if *open == 0 && b == c {
            return Some(index);
        }
    }
    None
}

fn erase_ranges(line: &mut Vec<u8>, limits: &[(usize, usize)]) {
    for &(first, last) in limits.iter().rev() {
        let end = (last + 1).min(line.len());
        if first < line.len() {
            line.drain(first..end.max(first));
        }
    }
    *line = truncate_spaces(line).to_vec();
}

/// `sStripNexusComments`: the comments of one line of a command (`[ ... ]`, nested ones inside
/// the outer pair), the open count carried over to the next line.
fn strip_comments(line: &mut Vec<u8>, open: &mut i32) {
    if is_blank(line) {
        return;
    }
    let mut limits: Vec<(usize, usize)> = Vec::new();
    let mut start = 0usize;
    let mut index = 0usize;
    while index < line.len() {
        let c = line[index];
        if c == b'[' {
            *open += 1;
            if *open == 1 {
                start = index;
            }
        } else if c == b']' {
            if *open == 1 {
                limits.push((start, index));
            }
            *open -= 1;
        }
        index += 1;
    }
    if *open != 0 {
        limits.push((start, index - 1));
    }
    erase_ranges(line, &limits);
}

/// `sStripCommentsOutsideCommand`: the comments before a command starts.
fn strip_comments_outside_command(line: &mut Vec<u8>, open: &mut i32, in_command: &mut bool) {
    if is_blank(line) {
        return;
    }
    let mut limits: Vec<(usize, usize)> = Vec::new();
    let mut start = 0usize;
    if !*in_command && *open == 0 && line[0] != b'[' {
        *in_command = true;
    }
    let len = line.len();
    for index in 0..len {
        let c = line[index];
        if *in_command {
            if c == b';' && *open == 0 {
                *in_command = false;
            }
            continue;
        }
        if c == b'[' {
            *open += 1;
            if *open == 1 {
                start = index;
            }
        } else if c == b']' {
            *open -= 1;
            if *open == 0 {
                limits.push((start, index));
            }
        } else if *open == 0 && !crate::ncbi_str::is_space(c) {
            *in_command = true;
        }
    }
    if *open != 0 && !*in_command {
        limits.push((start, len - 1));
    }
    erase_ranges(line, &limits);
}

/// `sStripNexusCommentsFromCommand`.
fn strip_comments_from_command(args: &mut Args) {
    let mut open = 0;
    let mut kept = Args::new();
    for mut line in std::mem::take(args) {
        strip_comments(&mut line.data, &mut open);
        if !line.data.is_empty() {
            kept.push(line);
        }
    }
    *args = kept;
}

struct Scanner {
    num_sequences: i32,
    sequence_size: i32,
    match_char: u8,
    missing_char: u8,
    gap_char: u8,
    in_block: bool,
    current_block: Vec<u8>,
    previous_sequin: bool,
    ids: Vec<Vec<u8>>,
    sequences: Vec<Vec<u8>>,
    deflines: Vec<Vec<u8>>,
}

fn fail<T>() -> Scan<T> {
    Err(ScanFailure::Failed)
}

fn undefined<T>() -> Scan<T> {
    Err(ScanFailure::Undefined(EMPTY_COMMAND))
}

impl Scanner {
    /// `xUnexpectedEndBlock`: a last line that is `end`. NCBI answers it with a throw in both of
    /// its branches (the error and the warning of the reporter), so the one value left is
    /// `false`; an empty command reads `back()` of an empty list.
    fn unexpected_end_block(&self, args: &Args) -> Scan<bool> {
        let Some(last) = args.last() else {
            return undefined();
        };
        let line = &last.data;
        let at = line.iter().rposition(|c| matches!(c, b' ' | b'\t'));
        let token = match at {
            None => &line[..],
            Some(at) => &line[at..],
        };
        if token.eq_ignore_ascii_case(b"end") {
            return fail();
        }
        Ok(false)
    }

    /// `xProcessCommand`.
    fn process_command(&mut self, tokens_: &Args) -> Scan<()> {
        let Some(front) = tokens_.first() else {
            return undefined();
        };
        let mut args: Args = tokens_.clone();
        let name: Vec<u8>;
        match find_first_of(&front.data, b" \t[", 0) {
            None => {
                name = front.data.clone();
                args.remove(0);
            }
            Some(end) => {
                name = front.data[..end].to_vec();
                args[0].data = truncate_spaces(&front.data[end..]).to_vec();
            }
        }
        let name_lower = name.to_ascii_lowercase();
        if name_lower == b"begin" {
            strip_comments_from_command(&mut args);
            let apply_end = self.unexpected_end_block(&args)?;
            self.begin_block(&args)?;
            if apply_end {
                self.end_block()?;
            }
            return Ok(());
        }
        if !self.in_block {
            return fail();
        }
        let block = self.current_block.to_ascii_lowercase();
        if block == b"ncbi" {
            return self.process_ncbi_block_command(&name_lower, &mut args);
        }
        if name_lower == b"end" {
            if !args.is_empty() {
                return fail();
            }
            return self.end_block();
        }
        if block == b"data" || block == b"characters" {
            return self.process_data_block_command(&name_lower, &mut args);
        }
        if block == b"taxa" {
            return self.process_taxa_block_command(&name_lower, &mut args);
        }
        Ok(())
    }

    fn process_data_block_command(&mut self, name: &[u8], args: &mut Args) -> Scan<()> {
        strip_comments_from_command(args);
        let unexpected_end = self.unexpected_end_block(args)?;
        match name {
            b"dimensions" => self.process_dimensions(args)?,
            b"format" => self.process_format(args),
            b"matrix" => self.process_matrix(args)?,
            _ => {}
        }
        if unexpected_end {
            self.end_block()?;
        }
        Ok(())
    }

    fn process_taxa_block_command(&mut self, name: &[u8], args: &mut Args) -> Scan<()> {
        strip_comments_from_command(args);
        let unexpected_end = self.unexpected_end_block(args)?;
        if name == b"dimensions" {
            self.process_dimensions(args)?;
        }
        if unexpected_end {
            self.end_block()?;
        }
        Ok(())
    }

    fn process_ncbi_block_command(&mut self, name: &[u8], args: &mut Args) -> Scan<()> {
        if name == b"end" {
            if !self.previous_sequin {
                return fail();
            }
            self.previous_sequin = false;
            return self.end_block();
        }
        let unexpected_end = self.unexpected_end_block(args)?;
        if name == b"sequin" {
            self.process_sequin(args)?;
            self.previous_sequin = true;
        } else {
            return fail();
        }
        if unexpected_end {
            self.previous_sequin = false;
            self.end_block()?;
        }
        Ok(())
    }

    fn process_sequin(&mut self, args: &Args) -> Scan<()> {
        for line in args {
            let (dummy, def) = process_defline(&line.data)?;
            if !dummy.is_empty() {
                return fail();
            }
            self.deflines.push(def);
        }
        Ok(())
    }

    fn begin_block(&mut self, args: &Args) -> Scan<()> {
        let Some(first) = args.first() else {
            return undefined();
        };
        if self.in_block {
            return fail();
        }
        self.in_block = true;
        self.current_block = first.data.clone();
        Ok(())
    }

    fn end_block(&mut self) -> Scan<()> {
        if !self.in_block {
            return fail();
        }
        self.in_block = false;
        self.current_block.clear();
        Ok(())
    }

    /// `xProcessMatrix`.
    fn process_matrix(&mut self, args: &Args) -> Scan<()> {
        let mut data_line_count = 0i32;
        let mut seq_count = 0i32;
        let mut block_line_length = 0i32;
        let mut sequence_char_count = 0i32;
        let mut max_seq_count = 0i32;
        for line in args {
            let toks = tokens(&line.data, b" \t");
            if toks.len() < 2 {
                return fail();
            }
            let seq_id = toks.at(0);
            if self.num_sequences == 0 && !self.ids.is_empty() && seq_id == self.ids[0].as_slice() {
                self.num_sequences = self.ids.len() as i32;
            }
            let in_first_block = self.num_sequences == 0 || data_line_count < self.num_sequences;
            seq_count = if in_first_block {
                data_line_count
            } else {
                data_line_count % self.num_sequences
            };
            max_seq_count = max_seq_count.max(seq_count);
            check_id(seq_id, &self.ids, seq_count as usize, in_first_block)?;
            if in_first_block {
                self.ids.push(seq_id.to_vec());
                self.sequences.push(Vec::new());
            }
            // The tokens after the id, joined: the line's data without its blanks.
            let rest: Vec<u8> = line_tokens_after_first(&line.data);
            let data_size = rest.len() as i32;
            if seq_count == 0 {
                sequence_char_count = sequence_char_count.wrapping_add(data_size);
                if sequence_char_count > self.sequence_size {
                    return fail();
                }
                block_line_length = data_size;
            } else if data_size != block_line_length {
                return fail();
            }
            self.sequences[seq_count as usize].extend_from_slice(&rest);
            data_line_count += 1;
        }
        if max_seq_count != self.num_sequences.wrapping_sub(1) {
            return fail();
        }
        if seq_count != max_seq_count {
            return fail();
        }
        if sequence_char_count != self.sequence_size {
            return fail();
        }
        Ok(())
    }

    /// `xGetKeyVal`: the value after `key` (possibly on a later line) and its line.
    fn key_value(args: &Args, key: &[u8]) -> (Vec<u8>, i32) {
        let mut key_pos: Option<usize> = None;
        let mut key_line = 0;
        let mut end_pos: Option<usize> = Some(0);
        for line in args {
            if key_pos.is_none() {
                key_pos = find_no_case(&line.data, key);
                if let Some(at) = key_pos {
                    key_line = line.num;
                    end_pos = find_first_of(&line.data, b" \t=", at);
                }
            }
            if key_pos.is_some() {
                if line.num != key_line {
                    end_pos = Some(0);
                }
                // `find_first_not_of(.., npos)` is npos.
                let val_pos = end_pos.and_then(|from| find_first_not_of(&line.data, b" \t=", from));
                if let Some(val_pos) = val_pos {
                    return match find_first_of(&line.data, b" \t\n;", val_pos) {
                        Some(end) => (line.data[val_pos..end].to_vec(), line.num),
                        None => (line.data[val_pos..].to_vec(), line.num),
                    };
                }
            }
        }
        (Vec::new(), -1)
    }

    /// `xProcessDimensions`.
    fn process_dimensions(&mut self, args: &Args) -> Scan<()> {
        if self.current_block.eq_ignore_ascii_case(b"characters") {
            // `xGetArgPos(args, "ntax")`: the first line that holds `ntax` (lower case).
            let found = args.iter().enumerate().find_map(|(at, line)| {
                let lower = line.data.to_ascii_lowercase();
                lower
                    .windows(4)
                    .position(|w| w == b"ntax")
                    .map(|pos| (at, pos))
            });
            if let Some((at, pos)) = found {
                let mut ntax_sub: Vec<u8> = Vec::new();
                let mut ntax_line_pos = pos;
                if ntax_line_pos == 0 && at != 0 {
                    ntax_sub = args[at - 1].data.clone();
                    ntax_line_pos += ntax_sub.len();
                }
                let mut found_error = true;
                const LIT: usize = 7;
                if ntax_line_pos > LIT {
                    ntax_sub.extend_from_slice(&args[at].data);
                    if let Some(end_prev) = find_last_not_of(&ntax_sub, b" \t", ntax_line_pos - 1) {
                        if end_prev >= LIT - 1
                            && ntax_sub[end_prev - (LIT - 1)..=end_prev]
                                .eq_ignore_ascii_case(b"newtaxa")
                        {
                            found_error = false;
                        }
                    }
                }
                if found_error {
                    return fail();
                }
            }
        }
        let (ntax, _) = Self::key_value(args, b"ntax");
        if !ntax.is_empty() {
            self.num_sequences = string_to_int(&ntax, 0).or_else(|_| fail())?;
        }
        let (nchar, _) = Self::key_value(args, b"nchar");
        if !nchar.is_empty() {
            self.sequence_size = string_to_int(&nchar, 0).or_else(|_| fail())?;
        }
        Ok(())
    }

    /// `xProcessFormat`: the first character of each value.
    fn process_format(&mut self, args: &Args) {
        let (missing, _) = Self::key_value(args, b"missing");
        let (gap, _) = Self::key_value(args, b"gap");
        let (matchchar, _) = Self::key_value(args, b"matchchar");
        if let Some(&c) = missing.first() {
            self.missing_char = c;
        }
        if let Some(&c) = gap.first() {
            self.gap_char = c;
        }
        if let Some(&c) = matchchar.first() {
            self.match_char = c;
        }
    }
}

/// The tokens of a matrix line after the first, joined (`Join(tokens.begin()+1, tokens.end(), "")`).
fn line_tokens_after_first(data: &[u8]) -> Vec<u8> {
    data.split(|c| matches!(c, b' ' | b'\t'))
        .filter(|token| !token.is_empty())
        .skip(1)
        .flatten()
        .copied()
        .collect()
}

/// `xImportAlignmentData`, `xAdjustSequenceInfo` and the base verification.
pub(crate) fn scan(input: &mut PeekAheadStream<'_>) -> Scan<Scanned> {
    let mut s = Scanner {
        num_sequences: 0,
        sequence_size: 0,
        match_char: 0,
        missing_char: 0,
        gap_char: 0,
        in_block: false,
        current_block: Vec::new(),
        previous_sequin: false,
        ids: Vec::new(),
        sequences: Vec::new(),
        deflines: Vec::new(),
    };
    let mut command_tokens: Args = Vec::new();
    let mut command_start = 0usize;
    let mut open = 0i32;
    let mut in_command = false;
    let mut first_token = true;
    while let Some((raw, line_count)) = input.read_line() {
        let mut line = truncate_spaces(raw).to_vec();
        if line.eq_ignore_ascii_case(b"#nexus") {
            if !first_token {
                return fail();
            }
            first_token = false;
            continue;
        }
        strip_comments_outside_command(&mut line, &mut open, &mut in_command);
        if line.is_empty() {
            continue;
        }
        let mut command_end = find_char_outside_comment(b';', &line, &mut open, 0);
        while let Some(end) = command_end {
            let substr = truncate_spaces(&line[command_start..end]).to_vec();
            if !substr.is_empty() {
                command_tokens.push(Line {
                    data: substr,
                    num: line_count,
                });
            }
            s.process_command(&command_tokens)?;
            command_tokens.clear();
            command_start = end + 1;
            command_end = find_char_outside_comment(b';', &line, &mut open, command_start);
        }
        if command_start < line.len() {
            command_tokens.push(Line {
                data: truncate_spaces(&line[command_start..]).to_vec(),
                num: line_count,
            });
        }
        command_start = 0;
    }
    if open > 0 {
        return fail();
    }
    if !command_tokens.is_empty() {
        return fail();
    }
    let mut info = SeqInfo::protein();
    if s.gap_char != 0 {
        info.gap = vec![s.gap_char];
    }
    if s.match_char != 0 {
        info.match_chars = vec![s.match_char];
    }
    if s.missing_char != 0 {
        info.missing = vec![s.missing_char];
    }
    verify_rows(&info, &s.sequences)?;
    Ok(Scanned {
        ids: s.ids,
        rows: s.sequences,
        deflines: s.deflines,
        info,
    })
}
